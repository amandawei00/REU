import csv

txt = "results0-30.txt"
cs = "results0-30.csv"

in_txt = csv.reader(open(txt, "rb"), delimiter='\t')
out_csv = csv.writer(open(cs, 'wb'))

out_csv.writerows(in_txt)
