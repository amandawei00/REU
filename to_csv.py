import csv

txt = "results_ATLAS.txt"
cs = "results_ATLAS.csv"

in_txt = csv.reader(open(txt, "rb"), delimiter='\t')
out_csv = csv.writer(open(cs, 'wb'))

out_csv.writerows(in_txt)
