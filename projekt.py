##Importy
import argparse
import httplib2
import pandas as pd
import json

##Parser argumentów z command line
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()  

parser = argparse.ArgumentParser(
                    prog='NAZWA PROGRAMU', ##IDK trzeba coś zdecydować
                    description='Program do zapytywanie MyVariants.info', ##Lepszy opis?
                    epilog='Bioinformatyka rok V',
                    formatter_class=SmartFormatter)
parser.add_argument('-i', '--input', default='', help='Sciezka do pliku VCF')
parser.add_argument('--id', default='', help='Wyszukuj warianty po ID; możliwe formaty:\n\t -RSid np. rs58991260 \n\t -konkretny SNP np: chr1:g.35367G>A \n\t -ENSBML gene ID np: ENSG00000113368')
parser.add_argument('-o', '--output', default='', help='Sciezka do zapisania raportu')
parser.add_argument("--show-na", action="store_true", help="Pokazuj warianty bez wpisów w bazach danych (domyslnie falsz)") #domyślnie fałsz, ale jak ktoś da w wywolaniu --test to włączy się prawda; do wykorzystania przy flag filtrowania
parser.add_argument("--rare", action="store_true", help="Pokazuje tylko rzadkie")
parser.add_argument("--pathogenic", action="store_true", help="Pokazuje tylko warianty klinicznie patogeniczne")
###Trzeba dodać tutaj inne flagi od funkcji filtrowania np
args = parser.parse_args()

#path = "testowy.vcf" ##Tymczasowa ścieżka do pliku, żeby nie trzeba było odpalać z CLI / do zakomentowania
#path = args.input

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
    for result in resultsJSON:
            if result.get("notfound", False):
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
                
                vcf = ""
                for key, value in result.get('vcf', {}).items():
                    vcf += f"{key}:{value}; "

                clinvar = ""
                for key, value in result.get('clinvar', {}).items():
                    if key == '_license':
                        continue
                    clinvar += f"{key}:{value}; "
                    

                snpeff = ""
                snpeff_type = type(result.get('snpeff',{}).get('ann', 'N/A'))
                if snpeff_type is dict:
                    for key, value in result['snpeff']['ann'].items():
                        snpeff += f"{key}:{value}; "
                        pass
                elif snpeff_type is list:
                    for annot in result.get('snpeff',{}).get('ann', 'N/A'):
                        for key, value in annot.items():
                            snpeff += f"{key}:{value}; "
                        snpeff += '|'
                else:
                    snpeff = result.get('snpeff', {}).get('ann', 'N/A')
                    
                dbsnp = ""
                rare = "N/A"
                for key, value in result.get('dbsnp', {}).items():
                    if key == '_license':
                        continue
                    if key == 'alleles':
                        rare = check_if_rare(value)
                    dbsnp += f"{key}:{value}; "

                if clinvar == "":
                    clinvar = "N/A"
                if dbsnp == "":
                    dbsnp = "N/A"
                if snpeff == "":
                    snpeff = "N/A"
                rows_to_append.append({'ID':id, 'SCORE':score, 'CHROM':chrom, 'START':start, 'END':end, 'OBSERVED':observed, 'VCF':vcf, 'CLINVAR':clinvar, 'SNPEFF':snpeff, 'DBSNP':dbsnp, 'RARE': rare})

    df = pd.DataFrame(rows_to_append)

    return df
    
##Zapisywanie raportu w HTML
#NATALIA, ZAPISUJESZ DO HTMLA W TYM MIEJSCU, CZUJ SIE WOLNA ZMIENIĆ WSZYSTKO CO CHCESZ
#https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_html.html
#mam nadzieje ze wam to tez dziala XD



def saveRaport(parsedJSON):
    if args.output == '':
        output = "raport.html"
    else:
        output = args.output + "\\raport.html"
    
    with open('skeleton.html', 'r') as skeleton:
        html = skeleton.read()

    table_html = parsedJSON.to_html(index=False)
    raport = html.replace("{table_html}", table_html)
    with open(output, "w") as file:
        file.write(raport)
        
    print(f"Report saved to {output}")
    return


def check_if_rare(alleles, bias=0.01):
    allele = alleles[-1]
    total_freq = 0
    db_counter = 0
    for db, freq in allele.get("freq", {}).items():
        if db != "":
            db_counter += 1
        total_freq += freq
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
    if args.rare:
        print("Showing only rare variants.")
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
        parsedJSON = parseJSON(resultsJSON)# <-- TUTAJ JEST DATAFRAME NATALIA, parseJSON(resultsJSON) ZWRACA DATAFRAME

        clinical_significance = []
        for result in resultsJSON:
            if result.get("notfound", False) and not args.show_na:
                continue 
            clinical_significance.append(ClinicalSignificance(result))

        parsedJSON["CLINICAL_SIGNIFICANCE"] = clinical_significance
        if args.rare:
            parsedJSON = parsedJSON[parsedJSON["RARE"]=="+"]

        if args.pathogenic:
            parsedJSON = parsedJSON[parsedJSON["CLINICAL_SIGNIFICANCE"].str.lower() == "pathogenic"]

        print("Saving raport...")
        saveRaport(parsedJSON) #Natalia modyfikuj plik skeleton.html. Tam, gdzie ma pojawiać się tabela wstaw <div id="wynik">{table_html}</div>
        
#Jak będzie trzeba coś jeszcze dodać po mojej stronie np handling flag jakiś czy coś wymyślicie to dajcie znać i ogarne ~Piotr