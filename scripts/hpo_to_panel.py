import requests, sys, json

#ensembl API: 
#https://grch37.rest.ensembl.org/
#https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/u1218/8_REST_API_Emily.pdf
#hgnc API:
#https://www.genenames.org/help/rest/#!/#tocAnchor-1-1-1

def fetch_endpoint(server, request, headers):
    r = requests.get(server+request, headers=headers)
    if not r.ok:
        return False
    if 'application/json' in headers.values():
        return r.json()
    else:
        return r.text

def search_ensmble(i):
    
    #grch37 url (set genome version in config later)
    server = "https://grch37.rest.ensembl.org"
    ext = "/lookup/id/{}?&expand=1".format(i)
    headers = {"Content-Type": con}
    response = fetch_endpoint(server, ext, headers)
    if response:
        #use only canonical chromosomes or return gene symbol
        if response["seq_region_name"] in canonical:
            return ",".join([ str(response[j]) for j in fields if j in response])
        else:
            try:
                not_found = response["display_name"]
            except:
                not_found = i
    #try to find gene name in hgnc
    else:
        response = search_hgnc(i,type="ensembl_gene_id")
        try:
            not_found =  response["response"]["docs"][0]["symbol"]
        except:
            not_found = i
    out =  "{},NotFound".format(not_found)
    return out

def search_hgnc(i,type="symbol"):
    server = "http://rest.genenames.org"
    ext = "/fetch/{}/{}".format(type,i)
    headers = {'Accept': con}
    response = fetch_endpoint(server, ext, headers)
    return response
                
def search_id(i):
    data = ''
    if i.startswith("ENS"):
        return search_ensmble(i)
    else:
        response = search_hgnc(i)
        if response:
            try:
                ensembl_gid = response['response']['docs'][0]['ensembl_gene_id']
                return search_ensmble(ensembl_gid)
            except:
                return "{},NotFound".format(i)
            

geneid = sys.argv[1]
con = "application/json"
fields = ["seq_region_name", "start", "end" ]
canonical =  list(map(str, range(1,23))) + ["X", "Y", "MT"]
print(search_id(geneid))
