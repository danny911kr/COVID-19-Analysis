import sys
import csv
import networkx as nx

# cathepsin (prefix)
# lectin-mannose-binding-1 (suffix)
# transmembrane-serine-protease-2 (TMPRSS2)
# furin (prefix)
# angiotensin-l-converting-enzyme-2 (ACE2)
target_list = []
seed_targets = {
    'ACE2': ['angiotensin I converting enzyme 2'],
    'furin': ['furin'],
    'TMPRSS2': ['Transmembrane protease serine 2', 'Transmembrane protease, serine 2', 'serine protease 102'],
    'cathepsin': ['cathepsin'],
    'mannose': ['lectin, mannose binding', 'mannose-binding lectin', 'mannose binding lectin', 'mannan-binding lectin', 'mannan binding lectin']
}

gene_dict = {}

target_gene_id = {
    'ACE2': [],
    'furin': [],
    'TMPRSS2': [],
    'cathepsin': [],
    'mannose': []
}

#from ctdbase
# coronavirus > enteritis, feline, gastroenteritis, covid-19, sars
disease_id = {
    # 'Coronavirus Infections': 'MESH:D003333',
    # 'Enteritis, Transmissible, of Turkeys': 'MESH:D004753',
    # 'Feline Infectious Peritonitis': 'MESH:D016766',
    # 'Gastroenteritis, Transmissible, of Swine': 'MESH:D005761',
    'COVID-19': 'MESH:C000657245',
    'SARS': 'MESH:D045169'
}
disease_id_invert = {
    # 'MESH:D003333':'Coronavirus Infections',
    # 'MESH:D004753':'Enteritis, Transmissible, of Turkeys',
    # 'MESH:D016766':'Feline Infectious Peritonitis',
    # 'MESH:D005761':'Gastroenteritis, Transmissible, of Swine',
    'MESH:C000657245':'COVID-19',
    'MESH:D045169':'SARS'
}
disease_id_to_graph_node = {
    # 'MESH:D003333':0,
    # 'MESH:D004753':1,
    # 'MESH:D016766':2,
    # 'MESH:D005761':3,
    'MESH:C000657245':0,
    'MESH:D045169':1
}

# # GeneSymbol, GeneName, GeneID, AltGeneIDs
with open('./KG/genes.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        for ace2 in seed_targets['ACE2']:
            if ace2.lower() in row['GeneName'].lower():
                target_gene_id['ACE2'].append(row['GeneID'])
                gene_dict[row['GeneID']] = row['GeneName']

        for furin in seed_targets['furin']:
            if furin.lower() in row['GeneName'].lower():
                target_gene_id['furin'].append(row['GeneID'])
                gene_dict[row['GeneID']] = row['GeneName']

        for tmprss2 in seed_targets['TMPRSS2']:
            if tmprss2.lower() in row['GeneName'].lower():
                target_gene_id['TMPRSS2'].append(row['GeneID'])
                gene_dict[row['GeneID']] = row['GeneName']

        for cathepsin in seed_targets['cathepsin']:
            if cathepsin.lower() in row['GeneName'].lower():
                target_gene_id['cathepsin'].append(row['GeneID'])
                gene_dict[row['GeneID']] = row['GeneName']

        for mannose in seed_targets['mannose']:
            if mannose.lower() in row['GeneName'].lower():
                target_gene_id['mannose'].append(row['GeneID'])
                gene_dict[row['GeneID']] = row['GeneName']

#
#frequency based intensity measurement
ace2 = {}
furin = {}
tmprss2 = {}
cathepsin = {}
mannose = {}

with open('./KG/genes_diseases_relation.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        for disease in disease_id:
            if disease_id[disease] in row['DiseaseID']:
                if row['GeneID'] in target_gene_id['ACE2']:
                    if row['GeneID'] in ace2:
                        if disease_id[disease] in ace2[row['GeneID']]:
                            ace2[row['GeneID']][disease_id[disease]].extend(row['pmids'].split('|'))
                        else:
                            ace2[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')
                    else:
                        ace2[row['GeneID']] = {}
                        ace2[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')

                if row['GeneID'] in target_gene_id['furin']:
                    if row['GeneID'] in furin:
                        if disease_id[disease] in furin[row['GeneID']]:
                            furin[row['GeneID']][disease_id[disease]].extend(row['pmids'].split('|'))
                        else:
                            furin[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')
                    else:
                        furin[row['GeneID']] = {}
                        furin[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')

                if row['GeneID'] in target_gene_id['TMPRSS2']:
                    if row['GeneID'] in tmprss2:
                        if disease_id[disease] in tmprss2[row['GeneID']]:
                            tmprss2[row['GeneID']][disease_id[disease]].extend(row['pmids'].split('|'))
                        else:
                            tmprss2[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')
                    else:
                        tmprss2[row['GeneID']] = {}
                        tmprss2[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')

                if row['GeneID'] in target_gene_id['cathepsin']:
                    if row['GeneID'] in cathepsin:
                        if disease_id[disease] in cathepsin[row['GeneID']]:
                            cathepsin[row['GeneID']][disease_id[disease]].extend(row['pmids'].split('|'))
                        else:
                            cathepsin[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')
                    else:
                        cathepsin[row['GeneID']] = {}
                        cathepsin[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')

                if row['GeneID'] in target_gene_id['mannose']:
                    if row['GeneID'] in mannose:
                        if disease_id[disease] in mannose[row['GeneID']]:
                            mannose[row['GeneID']][disease_id[disease]].extend(row['pmids'].split('|'))
                        else:
                            mannose[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')
                    else:
                        mannose[row['GeneID']] = {}
                        mannose[row['GeneID']][disease_id[disease]] = row['pmids'].split('|')


# ace2={'59272': {'MESH:C000657245': ['31996437', '32061198', '32081428', '32092392', '32117569', '32129518', '32132184', '32142651', '32149769', '32129518'], 'MESH:D045169': ['31996437', '32092392']}}
# furin={'23241': {'MESH:D045169': ['32125140']}}
# tmprss2={}
# cathepsin={'1508': {'MESH:C000657245': ['32020029', '32070753'], 'MESH:D045169': ['25666761']}, '1075': {'MESH:D045169': ['32125140']}, '1509': {'MESH:C000657245': ['32020029', '32070753'], 'MESH:D045169': ['15132791', '32125140']}, '1513': {'MESH:C000657245': ['32162456'], 'MESH:D045169': ['32125140']}, '1514': {'MESH:C000657245': ['32020029', '32070753', '32074550', '32133962', '32147628'], 'MESH:D045169': ['15351731']}, '321453': {'MESH:D045169': ['32125140']}, '554157': {'MESH:D045169': ['32125140']}}
# mannose={}

G = nx.Graph()

# G.add_node(0, label='Coronavirus Infections', modularity_class=0, sid='MESH:D003333')
# G.nodes[0]['viz'] = {'size': '140.0'}

# G.add_node(1, label='Enteritis, Transmissible, of Turkeys', modularity_class=1, sid='MESH:D004753')
# G.nodes[1]['viz'] = {'size': '100.0'}
# G.add_edge(0, 1)
#
# G.add_node(2, label='Feline Infectious Peritonitis', modularity_class=2, sid='MESH:D016766')
# G.nodes[2]['viz'] = {'size': '100.0'}
# G.add_edge(0, 2)
#
# G.add_node(3, label='Gastroenteritis, Transmissible, of Swine', modularity_class=3, sid='MESH:D005761')
# G.nodes[3]['viz'] = {'size': '100.0'}
# G.add_edge(0, 3)

G.add_node(0, label='SARS', modularity_class=0, sid='MESH:D045169')
G.nodes[0]['viz'] = {'size': '100.0'}
# G.add_edge(0, 1)

G.add_node(1, label='COVID-19', modularity_class=1, sid='MESH:C000657245')
G.nodes[1]['viz'] = {'size': '100.0'}
# G.add_edge(0, 1)

#https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=my_email@example.com&ids=32020029
if len(ace2) > 0:
    for a in ace2:
        G.add_node(len(G), label=gene_dict[a], modularity_class=2, sid=a)
        G.nodes[len(G) - 1]['viz'] = {'size': '40.0'}
        for ad in ace2[a]:
            G.add_edge(disease_id_to_graph_node[ad], len(G)-1, weight=len(ace2[a]))

if len(furin) > 0:
    for a in furin:
        G.add_node(len(G), label=gene_dict[a], modularity_class=3, sid=a)
        G.nodes[len(G)-1]['viz'] = {'size': '40.0'}
        for ad in furin[a]:
            G.add_edge(disease_id_to_graph_node[ad], len(G)-1, weight=len(furin[a]))

if len(tmprss2) > 0:
    for a in tmprss2:
        G.add_node(len(G), label=gene_dict[a], modularity_class=4, sid=a)
        G.nodes[len(G)-1]['viz'] = {'size': '40.0'}
        for ad in tmprss2[a]:
            G.add_edge(disease_id_to_graph_node[ad], len(G)-1, weight=len(tmprss2[a]))

if len(cathepsin) > 0:
    for a in cathepsin:
        G.add_node(len(G), label=gene_dict[a], modularity_class=5, sid=a)
        G.nodes[len(G)-1]['viz'] = {'size': '40.0'}
        for ad in cathepsin[a]:
            G.add_edge(disease_id_to_graph_node[ad], len(G) - 1, weight=len(cathepsin[a]))

if len(mannose) > 0:
    for a in mannose:
        G.add_node(len(G), label=gene_dict[a], modularity_class=6, sid=a)
        G.nodes[len(G)-1]['viz'] = {'size': '40.0'}
        for ad in mannose[a]:
            G.add_edge(disease_id_to_graph_node[ad], len(G) - 1, weight=len(mannose[a]))

nx.write_gexf(G, "corona_gene.gexf")

def invert_dict(dictfiles):
    inverted = {}
    for dictf in dictfiles:
        for gene_key in dictf:
            for disease_key in dictf[gene_key]:
                if disease_key in inverted:
                    if gene_key in inverted[disease_key]:
                        inverted[disease_key][gene_key].extend(dictf[gene_key][disease_key])
                    else:
                        inverted[disease_key][gene_key] = dictf[gene_key][disease_key]
                else:
                    inverted[disease_key] = {}
                    inverted[disease_key][gene_key] = dictf[gene_key][disease_key]
    return inverted

inverted = invert_dict([ace2,furin,tmprss2,cathepsin,mannose])
fileout = open("corona_gene_table.html", "w")

table = "<table class=\"tg\">\n"
table += "  <tr>\n"
table += "    <th class=\"tg-0lax\">{0}</th>\n".format('Disease')
table += "    <th class=\"tg-0lax\">{0}</th>\n".format('Gene')
table += "    <th class=\"tg-0lax\">{0}</th>\n".format('paper_link')
table += "  </tr>\n"

import requests
# Create the table's row data
for disease in inverted:
    merge_row = len(inverted[disease])
    for i, gene_key in enumerate(inverted[disease]):
        table += "  <tr class=\"tg-cly1\">\n"
        if i == 0:
            table += "    <td class=\"tg-cly1\" rowspan=\"{0}\">{1}</td>\n".format(merge_row, disease_id_invert[disease])
        table += "    <td class=\"tg-cly1\">{0}</td>\n".format(gene_dict[gene_key])
        table += "    <td class=\"tg-cly1\">"
        for pid in inverted[disease][gene_key]:
            r = requests.get('https://api.altmetric.com/v1/pmid/{0}'.format(pid))
            if r.status_code == 404:
                url = 'NONE'
            else:
                parsed = r.json()
                url = parsed['url']
                title = parsed['title']
            table += "<pre><a href=\"{0}\">{1}</a></pre>\n".format(url, title)
        table += "</td>\n"
        table += "  </tr>\n"
table += "</table>"

fileout.writelines(table)
fileout.close()
















