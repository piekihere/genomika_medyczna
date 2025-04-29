##Importy
import argparse
import httplib2
import pandas as pd
import json

##Parser argumentów z command line
parser = argparse.ArgumentParser(
                    prog='NAZWA PROGRAMU', ##IDK trzeba coś zdecydować
                    description='Program do zapytywanie MyVariants.info', ##Lepszy opis?
                    epilog='Bioinformatyka rok V')
parser.add_argument('-i', '--input', help='Sciezka do pliku VCF')
parser.add_argument('-o', '--output', default='', help='Sciezka do zapisania raportu')
parser.add_argument("--test", action="store_true", help="TESTOWY") #domyślnie fałsz, ale jak ktoś da w wywolaniu --test to włączy się prawda; do wykorzystania przy flag filtrowania
parser.add_argument("--minScore", type=int, defautl=0, help="TESTOWY") #Domyslnie 0; też do użycia przy filtorwaniu wynikow
###Trzeba dodać tutaj inne flagi od funkcji filtrowania np
args = parser.parse_args()

path = "testowy.vcf" ##Tymczasowa ścieżka do pliku, żeby nie trzeba było odpalać z CLI / do zakomentowania
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
                 'FORMAT':[],
                 'NORMAL':[],
                 'TUMOR':[]})
    
    try:
        with open(path, "r") as file:
            for line in file:
                if line[0] == "#":
                    continue
                values = line.split(sep="\t")
                values = [n.strip('\n') for n in values]
                vcfValues.loc[-1] = values
                vcfValues.index = vcfValues.index + 1
                vcfValues = vcfValues.sort_index()
            return vcfValues
        
    except FileNotFoundError:
        print("The file does not exist.")

##Tworzenie kwerendy dla MyVariants.info
#można dodać urozmaicanie kwerend może? https://docs.myvariant.info/en/latest/doc/variant_query_service.html#query-syntax
def makeQuery(vcf):
    ids = vcf['CHROM'] + ':g.' + vcf['POS'] + vcf['REF'] + '>' + vcf['ALT']
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
    
##Zapisywanie raportu w HTML
def saveRaport():
    output = args.output + "\\raport.html"
    print(output)
    return

##MAIN
if __name__=="__main__":
    vcf = readVCF(path)
    ###Można też przefiltorwać tutaj po paramaterach z samego pliku VCF i potem takiego dataframe'a przekazać do makeQuery https://docs.myvariant.info/en/latest/doc/variant_query_service.html#query-syntax
    # query = makeQuery(vcf)
    # results = askAPI(query)
    # if results == 1:
    #     print("MyVaraints.info server error")
    #     exit()
    # else:
    #     resultsJSON = json.loads(results)
    print(args.test)
    saveRaport()
#Jak będzie trzeba coś jeszcze dodać po mojej stronie np handling flag jakiś czy coś wymyślicie to dajcie znać i ogarne ~Piotr