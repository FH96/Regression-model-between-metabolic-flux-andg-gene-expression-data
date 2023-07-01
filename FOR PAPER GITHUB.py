import csv
import statistics
import cobra
import pandas
import sklearn.linear_model
from matplotlib import pyplot as plt
from sklearn.model_selection import LeaveOneOut

model=cobra.io.read_sbml_model("iJO1366.xml")

expression_values=pandas.read_csv("expression values.csv",index_col='genes')

flux=pandas.read_csv("Rnx-Flux relationship modified for MLR for each reaction.csv")

data=pandas.read_csv("D:/semesters/11/3/list of genes added-Rnx-Flux relationship modified for MLR for each reaction .csv")


rxns_f_list=pandas.read_csv("D:/semesters/11/first data ref 17/reactions.csv")

list_of_systematic_names_of_needed_genes=[]

for item in rxns_f_list:
    g_r_rule=model.reactions.get_by_id(item).gene_reaction_rule
    if 'or' in g_r_rule and 'and' in g_r_rule:
        or_splitted = g_r_rule.split(' or ')
        for i in or_splitted:
            if not 'and' in i:
                list_of_systematic_names_of_needed_genes.append(i)

            else:
                new_i = i.replace("( ", '')
                newer_i = new_i.replace(" )", '')

                list_of_systematic_names_of_needed_genes.extend(newer_i.split(' and '))

    elif g_r_rule=="":
        pass
        #print(item ,'has no gene and',g_r_rule)
    elif 'or' in g_r_rule:
        list_of_systematic_names_of_needed_genes.extend(g_r_rule.split(' or '))


    elif 'and' in g_r_rule:
       # print(item,'=',g_r_rule)
        list_of_systematic_names_of_needed_genes.extend(g_r_rule.split(' and '))


    else:

        list_of_systematic_names_of_needed_genes.append(g_r_rule)

'length= 126 '

#removing repetetive files:
list_of_systematic_names_of_needed_genes=list(dict.fromkeys
                                              (list_of_systematic_names_of_needed_genes))



list_of_systematic_names_of_needed_genes.to_csv("genes needed.csv")

column_of_reactions=[]

column_of_genes=[]
for n in range(flux.shape[0]):
    rxns=flux.iloc[n]["Rnx-Flux relationship"]
    rxns = rxns.replace("(", '')
    rxns = rxns.replace(")", '')
    rxns = rxns.replace(" ", '')

    list_of_rxns=[]
    if "+" in rxns:
        list_of_rxns=rxns.split("+")
    elif '-' in rxns:
        list_of_rxns=rxns.split("-")
    elif ',' in rxns:
        list_of_rxns=rxns.split(",")
    else:
        list_of_rxns.append(rxns)

    column_of_reactions.append(list_of_rxns)

    genes=[]
    for rxn in list_of_rxns:
        #print("list_of_rxns: ",list_of_rxns)
        #print("rxn: ",rxn)
        gene_rule=model.reactions.get_by_id(rxn).gene_reaction_rule

        if "or" in gene_rule:
            genes.extend(gene_rule.split(" or "))
        elif "and" in gene_rule:
            genes.extend(gene_rule.split(" and "))
        else:
            genes.append(gene_rule)
        genes=list(set(genes))

    column_of_genes.append(genes)

flux_and_genes=flux
flux_and_genes["list of reactions"]=column_of_reactions
flux_and_genes["list of genes"]=column_of_genes

#flux_and_genes.to_csv("list of genes added-Rnx-Flux relationship modified for MLR for each reaction .csv")
flux_and_genes.set_index("Flux Module Name (short)",inplace=True)

conditions=['Acetate','Fructose','Galactose','Glucose','Glycerol','Gluconate',
           'Pyruvate','Succinate']

pearson_for_reactions= {}
for reaction in flux_and_genes.index:

    Y = flux_and_genes.loc[reaction,conditions]
    X = expression_values.loc[flux_and_genes.loc[reaction,"list of genes"],conditions]
    X=X.T

    loocv = LeaveOneOut()

    y_test_for_pearson = []
    y_predicted_for_pearson = []

    for train_index, test_index in loocv.split(X):
        x_train = X.iloc[train_index]
        y_train = Y.iloc[train_index]

        x_test = X.iloc[test_index]
        y_test = Y.iloc[test_index]
        MLR_model = sklearn.linear_model.LinearRegression()


        MLR_model.fit(x_train, y_train)


        predicted = MLR_model.predict(x_test)
        y_test_for_pearson.append(y_test.iloc[0])
        y_predicted_for_pearson.append(predicted)

    y_test_for_pearson_df = pandas.DataFrame(y_test_for_pearson)
    y_predicted_for_pearson_df = pandas.DataFrame(y_predicted_for_pearson)

    #print('model score: ' , MLR_model.score(X,Y))
    #print("y_test_for_pearson_df: ",y_test_for_pearson_df)
    pearson = y_test_for_pearson_df.corrwith(y_predicted_for_pearson_df, axis=0)
    pearson_for_reactions[reaction] = pearson.iloc[0]


pearson_for_reactions_df=pandas.DataFrame(pearson_for_reactions,index=["pearson"])
pearson_for_reactions_df=pearson_for_reactions_df.T

#pearson_for_reactions_df.columns=["Flux Module Name (short)"]

#print("pearson_for_reactions as df: ",pearson_for_reactions_df)

