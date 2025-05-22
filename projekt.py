##Importy
import argparse
import httplib2
import pandas as pd
import json
import datetime

def wrap_collapsible(label: str, content: str) -> str:
    return f'''
    <details>
      <summary>{label}</summary>
      <div style="white-space:pre-wrap; padding: 5px;">{content}</div>
    </details>
    '''
#
#def auto_wrap_collapsible(field_name: str, value: str, limit: int = 80) -> str:
#    if isinstance(value, str) and len(value) > limit:
#        return wrap_collapsible(f"Show {field_name}", value)
#    return value
def prettify_nested(value, indent=0):
    if isinstance(value, dict):
        html = ""
        for k, v in value.items():
            html += f"<details><summary>{k}</summary>{prettify_nested(v, indent + 1)}</details>"
        return html
    elif isinstance(value, list):
        html = ""
        for item in value:
            html += f"<div style='margin-left:{indent * 15}px'>{prettify_nested(item, indent + 1)}</div>"
        return html
    else:
        cleaned = str(value).replace("\\n", "<br>").replace("\n", "<br>").replace("<br><br>", "<br>")

        return f"<div style='margin-left:{indent * 15}px'>{cleaned}</div>"


def auto_wrap_collapsible(field_name: str, value: str, limit: int = 80) -> str:
    try:
        parsed = json.loads(value)
        if isinstance(parsed, (dict, list)):
            content = prettify_nested(parsed)
            return wrap_collapsible(f"Show {field_name}", content)
    except Exception:
        pass

    clean_value = str(value).replace("\\n", "<br>").replace("\n", "<br>").replace("<br><br>", "<br>")
    if isinstance(value, str) and len(value) > limit:
        return wrap_collapsible(f"Show {field_name}", clean_value)
    return clean_value



def df_to_html(df):
    df_copy = df.copy()
    for col in ["VCF", "CLINVAR", "SNPEFF", "DBSNP"]:
        if col in df_copy.columns:
            df_copy[col] = df_copy[col].apply(lambda x: auto_wrap_collapsible(col, x))
    return df_copy.to_html(classes="display", escape=False, index=False)

##Parser argumentów z command line
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()  

parser = argparse.ArgumentParser(
                    prog='Varianinator', ##IDK trzeba coś zdecydować
                    description='Search ClinVar, DBSNP and SnpEff for described variant in VCF file or by variant ID.',
                    epilog='Bioinformatyka rok V',
                    formatter_class=SmartFormatter)
parser.add_argument('-i', '--input', default='', help='Path to VCF file')
parser.add_argument('--id', default='', help='Search for variants by ID. Possible formats:\n\t -RSid e.g. rs58991260 \n\t -specific SNP e.g. chr1:g.35367G>A \n\t -ENSBML gene ID e.g. ENSG00000113368')
parser.add_argument('-o', '--output', default='', help='Path to output file.')
parser.add_argument("--show-na", action="store_true", help="Show variants even if there is no information for them in the databases") #domyślnie fałsz, ale jak ktoś da w wywolaniu --test to włączy się prawda; do wykorzystania przy flag filtrowania
parser.add_argument("--rare", action="store_true", help="Save only rare variants")
parser.add_argument("--pathogenic", action="store_true", help="Save only pathogenic variants.")
###Trzeba dodać tutaj inne flagi od funkcji filtrowania np
args = parser.parse_args()

##Wczytywanie pliku
def readVCF(path):
    vcfValues = pd.DataFrame({'CHROM':[], 
                 'POS':[],
                 'ID':[],
                 'REF':[],
                 'ALT':[],
                 'QUAL':[],
                 'FILTER':[],
                 'INFO':[],
                 'FORMAT':[]
                 })
    
    try:
        with open(path, "r") as file:
            for line in file:
                if line[0] == "#":
                    continue
                values = line.split(sep="\t")
                values = [n.strip('\n') for n in values[:9]]
                vcfValues.loc[-1] = values
                vcfValues.index = vcfValues.index + 1
                vcfValues = vcfValues.sort_index()
            return vcfValues
        
    except FileNotFoundError:
        print("The file does not exist.")
        exit()

##Tworzenie kwerendy dla MyVariants.info
def makeQuery(values, byID):
    if byID:
        if values[:2] == "rs":
            query = 'q='+ values + '&scopes=dbsnp.rsid'
        elif values[:4] == "ENSG":
            query = 'q=' + values + '&scopes=dbnsfp.ensembl.geneid'
        elif ">" in values:
            query = 'q=' + values
        else:
            print("Invalid ID syntax. Please use one of the following:\n\t -RSid np. rs58991260 \n\t -konkretny SNP np: chr1:g.35367G>A \n\t -ENSBML gene ID np: ENSG00000113368")
            #parser.print_help()
            exit()
    else:
        ids = values['CHROM'] + ':g.' + values['POS'] + values['REF'] + '>' + values['ALT']
        query = 'q=' + ','.join(ids)
    return query

##Wysyłanie i odbieranie zapytania
def askAPI(query):
    h = httplib2.Http()
    headers = {'content-type': 'application/x-www-form-urlencoded'}
    try:
        response, content = h.request('http://myvariant.info/v1/query', 'POST', query, headers=headers)
    except httplib2.ServerNotFoundError:
        print("No internet connection")
        return 1
    if response.status != 200:
        return 1
    else:
        contentDecoded = content.decode()
        return contentDecoded

#Parsowanie JSONa na dataframe
def parseJSON(resultsJSON):
    rows_to_append = []
    na_num = 0
    for result in resultsJSON:
            if result.get("notfound", False):
                na_num += 1
                if args.show_na == True:
                    id = result.get('query', 'N/A')
                    rows_to_append.append({'ID':id, 'SCORE':'N/A', 'CHROM':'N/A', 'START':'N/A', 'END':'N/A', 'OBSERVED':'N/A', 'VCF':'N/A', 'CLINVAR':'N/A', 'SNPEFF':'N/A', 'DBSNP':'N/A', 'RARE':'N/A'})
                else:
                    continue
            else:
                id = result.get('_id', 'N/A')
                score = result.get('_score', 'N/A')
                chrom = result.get('chrom', 'N/A')
                start = result.get('hg19', {}).get('start', 'N/A')
                end   = result.get('hg19', {}).get('end', 'N/A')
                observed = result.get('observed', 'N/A')
                
                vcf_dict = result.get('vcf', {})
                vcf = json.dumps(vcf_dict)

                clinvar_dict = {k: v for k, v in result.get('clinvar', {}).items() if k != '_license'}
                clinvar = json.dumps(clinvar_dict)
                if not bool(clinvar_dict):
                    clinvar = "N/A"

                snpeff_data = result.get('snpeff', {}).get('ann', 'N/A')
                snpeff = json.dumps(snpeff_data)
                if not bool(snpeff_data):
                    clinvar = "N/A"

                dbsnp_dict = {k: v for k, v in result.get('dbsnp', {}).items() if k != '_license'}
                if 'alleles' in dbsnp_dict:
                    rare = check_if_rare(dbsnp_dict['alleles'])
                else:
                    rare = "N/A"
                dbsnp = json.dumps(dbsnp_dict)
                if not bool(dbsnp_dict):
                    dbsnp = "N/A"


                if clinvar == "":
                    clinvar = "N/A"
                if dbsnp == "":
                    dbsnp = "N/A"
                if snpeff == "":
                    snpeff = "N/A"
                rows_to_append.append({'ID':id, 'SCORE':score, 'CHROM':chrom, 'START':start, 'END':end, 'OBSERVED':observed, 'VCF':vcf, 'CLINVAR':clinvar, 'SNPEFF':snpeff, 'DBSNP':dbsnp, 'RARE': rare})

    df = pd.DataFrame(rows_to_append)

    return df, na_num
    
##Zapisywanie raportu w HTML
def saveReport(parsedJSON):
    if args.output == '':
        output = "report.html"
    else:
        if args.output[-5:] == ".html":
            output = args.output
        else:
             output = args.output + ".html"
    
    with open('skeleton.html', 'r') as skeleton:
        html = skeleton.read()
    
    time = str(datetime.datetime.now())[:-10]
    table_html = df_to_html(parsedJSON)
    html = html.replace("{table_html}", table_html)
    html = html.replace("{time}", time)

    with open(output, "w") as file:
        file.write(html)
        
    print(f"Report saved to {output}")
    return

##SPrawdzanie rzadkości
def check_if_rare(alleles, bias=0.01):
    allele = alleles[-1]
    total_freq = 0
    db_counter = 0
    for db, freq in allele.get("freq", {}).items():
        if db != "":
            db_counter += 1
        total_freq += freq
    if db_counter != 0:
        avg = total_freq / db_counter
        if avg > bias:
            return "-"
        else:
            return "+"
    return "N/A"

# CLING SIG
def ClinicalSignificance(result):
    return result.get("clinvar", {}).get("rcv", {}).get("clinical_significance", "N/A")

##MAIN
if __name__=="__main__":
    if args.show_na:
        print("Showing empty results enabled.")
    if args.input == "" and args.id == "":
        print("Please use --input for VCF file or --id for searching variants by id!")
        exit()
    if args.input != "":
        print("Reading VCF file...")
        byID = False
        values = readVCF(args.input)
    elif args.id != "":
        print("Searching for variants by ID...")
        byID = True
        values = args.id

    print("Generating queries...")
    query = makeQuery(values, byID)
    print("Connecting to server...")
    results = askAPI(query)
    if results == 1:
        print("MyVaraints.info server error")
        exit()
    else:
        print("Fetching results...")
        resultsJSON = json.loads(results)
        if byID:
            if 'notfound' in resultsJSON[0].keys():
                print("Cannot find variant with given ID.")
                exit()
        parsedJSON, na_num = parseJSON(resultsJSON)# <-- TUTAJ JEST DATAFRAME NATALIA, parseJSON(resultsJSON) ZWRACA DATAFRAME
        print(f'Data for {na_num} varianst was not available.')

        clinical_significance = []
        for result in resultsJSON:
            if result.get("notfound", False) and not args.show_na:
                continue 
            clinical_significance.append(ClinicalSignificance(result))
        parsedJSON["CLINICAL_SIGNIFICANCE"] = clinical_significance

        if args.rare:
            before = len(parsedJSON)
            parsedJSON = parsedJSON[parsedJSON["RARE"]=="+"]
            after = len(parsedJSON)
            if parsedJSON.empty:
                print("No rare variants found!")
                exit()
            else:
                print(f'From {before}, {after} were considerd rare. {before - after} entries were omitted.')
                print("Saving only rare variants.")

        if args.pathogenic:
            before = len(parsedJSON)
            parsedJSON = parsedJSON[parsedJSON["CLINICAL_SIGNIFICANCE"].str.lower().isin(["pathogenic", "likely pathogenic"])]
            after = len(parsedJSON)
            if parsedJSON.empty:
                print("No pathogenic variants found!")
                exit()
            else:
                print(f'From {before}, {after} were considerd pathogenic. {before - after} entries were omitted.')
                print("Saving only pathogenic variants.")

        print("Saving report...")
        saveReport(parsedJSON) #Natalia modyfikuj plik skeleton.html. Tam, gdzie ma pojawiać się tabela wstaw <div id="wynik">{table_html}</div>
        
#Jak będzie trzeba coś jeszcze dodać po mojej stronie np handling flag jakiś czy coś wymyślicie to dajcie znać i ogarne ~Piotr
