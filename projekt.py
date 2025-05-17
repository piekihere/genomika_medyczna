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
parser.add_argument("--minScore", type=int, default=0, help="TESTOWY") #Domyslnie 0; też do użycia przy filtorwaniu wynikow
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

#Parsowanie JSONa na dataframe
def parseJSON(resultsJSON):
    rows_to_append = []
    for result in resultsJSON:
            if result.get("notfound", False):
                id = result.get('query', 'N/A')
                rows_to_append.append({'ID':id, 'SCORE':'N/A', 'CHROM':'N/A', 'START':'N/A', 'END':'N/A', 'OBSERVED':'N/A', 'VCF':'N/A', 'SNPEFF':'N/A', 'DBSNP':'N/A'})
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
                for key, value in result.get('dbsnp', {}).items():
                    if key == '_license':
                        continue
                    dbsnp += f"{key}:{value}; "
                    
                rows_to_append.append({'ID':id, 'SCORE':score, 'CHROM':chrom, 'START':start, 'END':end, 'OBSERVED':observed, 'VCF':vcf, 'SNPEFF':snpeff, 'DBSNP':dbsnp})


    df = pd.DataFrame(rows_to_append)

    return df
    
##Zapisywanie raportu w HTML
#NATALIA, ZAPISUJESZ DO HTMLA W TYM MIEJSCU, CZUJ SIE WOLNA ZMIENIĆ WSZYSTKO CO CHCESZ
#https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_html.html
#mam nadzieje ze wam to tez dziala XD
    table_html = df.to_html(index=False, border=1, classes='display')
    full_html = f"""
<!DOCTYPE html>
<html lang="pl">
<head>
    <meta charset="UTF-8">
    <title>Projekt</title>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <style>
        body {{
            text-align: center;
            font-family: Comic Sans MS;
            background-image: url('https://media1.giphy.com/media/v1.Y2lkPTc5MGI3NjExMzZ0czE4eW4wZWsxNWhiaGh3cjhoZHFxOW9leWd5eDZ2cGZvbnczdCZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/ECWLeGxcKasCc/giphy.gif');
            background-size: cover;
            color: #000;
        }}
        table.display {{
            background-color: rgba(255,255,255,0.8);
            margin: 20px auto;
            border-collapse: collapse;
        }}
        th, td {{ padding: 10px; border: 1px solid #aaa; }}
    </style>
    <script src="https://code.jquery.com/jquery-3.7.0.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
</head>
<body>
    <h1>PROJEKTTTT</h1>
    <img src="https://media.giphy.com/media/13borq7Zo2kulO/giphy.gif" width="300" alt="Jednorożec">
    <img src="https://media1.giphy.com/media/v1.Y2lkPTc5MGI3NjExY245NG84OXE5eDk1ZW5mOG1xM3E5bHR6emR6OHRjNmtxNHViajdhZyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/2A75RyXVzzSI2bx4Gj/giphy.gif" width="300" alt="Tęcza">

    <div>
        {table_html}
    </div>

    <script>
        $(document).ready(function() {{
            $('table.display').DataTable({{
                paging: true,
                searching: true,
                ordering: true,
                info: true
            }});
        }});
    </script>

    <audio autoplay loop>
        <source src="magiczna_muzyka.mp3" type="audio/mpeg">
    </audio>

    <script>
        function spawnRandomImage() {{
            const img = document.createElement('img');
            img.src = 'phd.png';
            img.style.position = 'fixed';
            img.style.width = '100px'; img.style.opacity = '0';
            const x = Math.random()*(window.innerWidth-100);
            const y = Math.random()*(window.innerHeight-100);
            img.style.left = x+'px'; img.style.top = y+'px';
            document.body.appendChild(img);
            setTimeout(() => img.style.opacity='1', 100);
            setTimeout(() => img.style.opacity='0', 3000);
            setTimeout(() => img.remove(), 3000);
        }}
        setInterval(spawnRandomImage, 3000);
    </script>

    <script>
        const mucha = document.createElement("img");
        mucha.src = "mucha.png";
        mucha.style.position = "fixed";
        mucha.style.width = "80px";
        mucha.style.zIndex = 9999;
        mucha.style.top = "50%";
        mucha.style.left = "50%";
        mucha.style.transition = "top 1s linear, left 1s linear";
        document.body.appendChild(mucha);

        function moveFly() {{
            const maxX = window.innerWidth - 80;
            const maxY = window.innerHeight - 80;
            const x = Math.random() * maxX;
            const y = Math.random() * maxY;
            mucha.style.left = x + "px";
            mucha.style.top = y + "px";
        }}

        setInterval(moveFly, 1000);
    </script>
</body>
</html>
"""

def saveRaport(parsedJSON):
    if args.output == '':
        output = "raport.html"
    else:
        output = args.output + "\\raport.html"
        
    with open(output, "w") as file:
        file.write(parsedJSON.to_html(index=False))
        
    print(f"Report saved to {output}")
    return

##MAIN
if __name__=="__main__":
    vcf = readVCF(path)
    ###Można też przefiltorwać tutaj po paramaterach z samego pliku VCF i potem takiego dataframe'a przekazać do makeQuery https://docs.myvariant.info/en/latest/doc/variant_query_service.html#query-syntax
    query = makeQuery(vcf)
    results = askAPI(query)
    if results == 1:
        print("MyVaraints.info server error")
        exit()
    else:
        resultsJSON = json.loads(results)
        parsedJSON = parseJSON(resultsJSON)# <-- TUTAJ JEST DATAFRAME NATALIA, parseJSON(resultsJSON) ZWRACA DATAFRAME
        saveRaport(parsedJSON) # <-- TUTAJ JEST ZAPISYWANIE DO HTML NATALIA, saveRaport(parsedJSON) ZAPISUJE DO HTML
        
#Jak będzie trzeba coś jeszcze dodać po mojej stronie np handling flag jakiś czy coś wymyślicie to dajcie znać i ogarne ~Piotr
