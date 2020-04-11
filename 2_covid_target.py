import sys
import csv
import networkx as nx

def parse_arguments(parser):
    parser.add_argument('--disease', type=str, default="covid-19")
    parser.add_argument('--chemical', type=str, default="")

    args = parser.parse_args()
    for k in args.__dict__:
        print(k + ": " + str(args.__dict__[k]))
    return args


# cathepsin (prefix)
# lectin-mannose-binding-1 (suffix)
# transmembrane-serine-protease-2 (TMPRSS2)
# furin (prefix)
# angiotensin-l-converting-enzyme-2 (ACE2)
target_list = []
seed_targets = {
    'ACE2': ['angiotensin I converting enzyme 2'], # include
    'furin': ['furin'], # lower, startswith
    'TMPRSS2': ['Transmembrane protease serine 2', 'Transmembrane protease, serine 2', 'serine protease 102'], # lower, exact
    'cathepsin': ['cathepsin'], #lower, startswith
    'mannose': ['lectin, mannose binding', 'mannose-binding lectin', 'mannose binding lectin', 'mannan-binding lectin', 'mannan binding lectin']
}

gene_dict = {}

# target_gene_id = {
#     'ACE2': [59272, 108707435, 100037182],
#     'furin': [47220, 32604, 5045, 566557, 566591, 398345, 397747, 24595193, 24595194, 10905857],
#     'TMPRSS2': [36342198, 100195190, 10907058, 100121542],
#     'cathepsin': [100195493, 100144796, 100144782, 100144784, 100144785, 100144791, 100144778, 100144794, 100144781, 100144790, 100144775, 100144795, 100144788, 100144787, 100144780, 692768, 100136231, 100528180, 100528797, 100286607, 100196462, 180111, 171829, 567046, 117066, 58518, 56092, 56094, 290972, 5476, 494810, 1508, 32341, 406645, 569298, 380102, 379257, 1075, 380203, 108710023, 1509, 777956, 443721, 1510, 417848, 373572, 373573, 8722, 108713957, 1511, 1512, 100036949, 26898, 1513, 100162831, 100127265, 1514, 436641, 484140, 321453, 30443, 449826, 444163, 70202, 64139, 1519, 108699249, 104002, 408201, 56835, 1520, 393398, 496647, 554157, 337572, 266984, 380516, 446505, 1515, 447313, 1521, 108705476, 1522, 432187, 494800, 36335937, 36336802, 36338414, 36343353, 36343838, 9949849, 24589432, 24594948, 24595674, 24595728, 24596741, 8342937, 8341971, 8347035, 36373312, 36373313, 36373450, 36375905, 36376055, 36377914, 36378791, 36378813, 36379415, 36379827, 36381286, 36381288, 36381290, 36383250, 36383654, 36383981, 10911327, 10908163, 10908162, 10908802, 10903649, 10905796, 10902555, 10912974],
#     'mannose': [3998, 79748, 108718427, 394295, 10960, 81562, 17194, 50639, 8512, 4153, 398413, 108696362, 373654, 5648, 10747]
# }
target_gene_id = {
    'ACE2': [],
    'furin': [],
    'TMPRSS2': [],
    'cathepsin': [],
    'mannose': []
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



# with open('./KG/diseases.csv', newline='') as csvfile:
#     reader = csv.DictReader(csvfile, delimiter='\t')
#     for row in reader:
#         if 'covid' in row['DiseaseName'].lower():
#             print(row['DiseaseName'], row['DiseaseID'])
disease_id = {
    'COVID-19': ['MESH:C000657245']
}

#frequency based intensity measurement
ace2 = {}
furin = {}
tmprss2 = {}
cathepsin = {}
mannose = {}

total = 0

with open('./KG/genes_diseases_relation.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        if disease_id['COVID-19'][0] in row['DiseaseID']:
            if row['GeneID'] in target_gene_id['ACE2']:
                if row['GeneID'] in ace2:
                    ace2[row['GeneID']] += 1
                else:
                    ace2[row['GeneID']] = 1
                total += 1

            if row['GeneID'] in target_gene_id['furin']:
                if row['GeneID'] in furin:
                    furin[row['GeneID']] += 1
                else:
                    furin[row['GeneID']] = 1
                total += 1
            if row['GeneID'] in target_gene_id['TMPRSS2']:
                if row['GeneID'] in tmprss2:
                    tmprss2[row['GeneID']] += 1
                else:
                    tmprss2[row['GeneID']] = 1
                total += 1
            if row['GeneID'] in target_gene_id['cathepsin']:
                if row['GeneID'] in cathepsin:
                    cathepsin[row['GeneID']] += 1
                else:
                    cathepsin[row['GeneID']] = 1
                total += 1
            if row['GeneID'] in target_gene_id['mannose']:
                if row['GeneID'] in mannose:
                    mannose[row['GeneID']] += 1
                else:
                    mannose[row['GeneID']] = 1
                total += 1

print(ace2)
print(furin)
print(tmprss2)
print(cathepsin)
print(mannose)

G = nx.Graph()

G.add_node(0, label='ACE2', modularity_class='0')
G.add_node(1, label='furin', modularity_class='1')
G.add_node(2, label='TMPRSS2', modularity_class='2')
G.add_node(3, label='cathepsin', modularity_class='3')
G.add_node(4, label='mannose', modularity_class='4')
G.add_node(5, label='COVID-19', modularity_class='5')

for a in ace2:
    ace2[a]




nx.write_gexf(G, "test.gexf")
