__author__ = 'emil'
import os
import zipfile
from six import BytesIO
import re
import bisect
from collections import defaultdict
import csv
import operator
import numpy as np
import pandas as pd

ENDMARGIN = 5

def find_le(a, x, **hi):
    'Find rightmost value less than or equal to x'
    i = bisect.bisect_right(a, x, **hi)
    if i:
        return i-1
    raise ValueError


def find_ge(a, x, **lo):
    'Find leftmost item greater than or equal to x'
    i = bisect.bisect_left(a, x, **lo)
    if i != len(a):
        return i
    raise ValueError


def kit_generator(maxn):
    with zipfile.ZipFile('Archive.zip') as zf:
        fl = zf.filelist
        yielded = 0
        for member in fl:
            if 'MACOSX' in member.filename:
                continue
            i = member.filename.split('_')[0] #Was int($)
            if maxn is not None and yielded > maxn:
                return

            with zf.open(member) as fp:

                zf_b = BytesIO(fp.read())
                zf_b.seek(0)
                sub_zf = zipfile.ZipFile(zf_b)
                try:
                    variant_fp = sub_zf.open('variants.vcf')
                except KeyError:
                    print("No variants.vcf in kit {0}".format(i))
                    
                regions_fp = sub_zf.open('regions.bed')
                yield variant_fp, regions_fp, i
                yielded += 1
                #raise StopIteration
                variant_fp.close()
                regions_fp.close()

snp_set = set()
snp_per_kit = defaultdict(set)
regions_per_kit = dict()

snp_ex_regx = re.compile(b'(\d+)\t\.\t([ACGT])\t([AGCT])\t\d+\.\d+\tPASS')

for variant, regions, i in kit_generator(None):
    kit_snp = snp_per_kit[i]
    try:
        var_row = variant.readline()
    except zipfile.BadZipFile:
        print("Skipping kit {0}".format(i))
        continue

    while b'chrY\t' not in var_row:
        var_row = variant.readline()

    while b'chrY\t' in var_row:
        m = snp_ex_regx.search(var_row)
        if m:
            if m.group(2) != m.group(3):
                snp_set.add(int(m.group(1)))
                kit_snp.add(int(m.group(1)))
                
                
        try:
            var_row = variant.readline()
        except zipfile.BadZipFile:
            print("Bad read in kit {0} line {1}".format(i,var_row))
            continue

    a,b = tuple(zip(*tuple(line.split(b'\t')[1:] for line in regions)))
    start = [int(_a) for _a in a]
    stop = [int(_b[:-1]) for _b in b]
    regions_per_kit[i] = (start, stop)

snp_list = sorted(int(snp) for snp in snp_set)
out = list()
no_of_positives = defaultdict(int)
n_snp = len(snp_list)
n_kits = len(regions_per_kit)
M = np.zeros((n_snp, n_kits))


for i, regs in regions_per_kit.items():
    snp_seen_in_kit = set()
    last_lo = 0
    for start, stop in zip(*regs):
        try:
            idx = find_ge(snp_list, start, lo=last_lo), find_le(snp_list, stop) + 1
            last_lo = idx[0]
            snp_in_region = snp_list[idx[0]:idx[1]]
        except ValueError:
            continue
        kit_snp_list = snp_per_kit[i]
        for snp in snp_in_region:
            snp_seen_in_kit.add(snp)
            if snp in kit_snp_list:
                out.append((i, snp, 'pos'))
                no_of_positives[snp] += 1
            elif snp - start > ENDMARGIN and stop - snp > ENDMARGIN:
                out.append((i, snp, 'neg'))
            else:
                out.append((i, snp, 'cb'))
    not_seen = snp_set.difference(snp_seen_in_kit)
    out.extend([(i, snp, 'nc') for snp in not_seen])

out.sort(key=lambda v: v[1])
out.sort(key=lambda v: v[0])
out_copy = list(out)
out_copy.reverse()

enumeration = {'pos': 1, 'neg': -1, 'cb': 0, 'nc': np.nan}
snp_map = list()
kit_map = list()

for i in range(n_kits):
    for j in range(n_snp):
        _out = out_copy.pop()
        #print("M[{0},{1}] = {2}".format(j,i, _out))
        M[j, i] = enumeration[_out[-1]]
        if i == 0:
            snp_map.append(_out[1])
    kit_map.append(_out[0])
#print(M)

np.save('allSNPs.npy',M)
## Creating output

df = pd.DataFrame(M, index=snp_map, columns=kit_map)
df.to_pickle('allSNPs.pkl')

with open('allSNPs.csv', 'w') as fp:
    writer = csv.writer(fp)
    writer.writerows(out)
    
with open('frequency.csv', 'w') as fp:
    writer = csv.writer(fp)
    writer.writerows(no_of_positives.items())
    
sorted_no_of_positives = sorted(no_of_positives.items(), key=operator.itemgetter(1), reverse=True)
with open('frequencysorted.csv', 'w') as fp:
    writer = csv.writer(fp)
    writer.writerows(sorted_no_of_positives)
