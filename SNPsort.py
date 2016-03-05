__author__ = 'pjp'
import os
import csv
from collections import defaultdict
import operator
import pandas as pd
import numpy as np
import cProfile
import sys


class SNPPrototype:
    def __init__(self,position):
        self.position = position
        self.primaryName = str(position)
        self.names = []
    def addName(self,name):
        if name not in self.names:
            self.names.append(name)
            if self.primaryName == str(self.position):
                self.primaryName = name
    def setPrimaryName(self,name):
        if name not in self.names:
            self.names.append(name)
        self.primaryName = name
    def findNames(self,SNPnamesAll):
        if str(self.position) in SNPnamesAll:
            for name in SNPnamesAll[str(self.position)]:
                self.addName(name)
                    
                    
class SNPequivalents:
    def __init__(self):
        self.primarySNP = 0
        self.name = ""
        self.SNPs = []
        self.posgroup = []
        self.plusgroup = [] 
        self.subclades = []
        self.neggroup = []
        self.foundAllSubclades = False
        self.failedQualityCheck = False
    def setPrimarySNP(self,primarySNP):
        self.primarySNP = primarySNP
        self.name = primarySNP.primaryName    
    def addSNP(self,SNP):
        if self.primarySNP == 0:
            self.setPrimarySNP(SNP)
        if SNP.primaryName != str(SNP.position):
            if self.name == str(self.primarySNP.position):
                self.setPrimarySNP(SNP)
            else:
                if SNP.position > 8000000 and SNP.position < self.primarySNP.position:
                    self.setPrimarySNP(SNP)
        self.SNPs.append(SNP)
    def addPosGroup(self,posgroup):
        self.posgroup = posgroup
    def addPlusGroup(self,plusgroup):
        self.plusgroup = plusgroup
    def addNegGroup(self,neggroup):
        self.neggroup = neggroup
    def addSubclade(self,subclade):
        self.subclades.append(subclade)
        
    def findSubcladeCandidate(self, df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll, skipSNPs, maxcount):
        
        ingroup = []
        inposplusgroup = []
        if len(self.subclades) == 0:
            ingroup = self.posgroup
            inposplusgroup = list(self.posgroup)         
        else:
            ingroup = self.subclades[len(self.subclades)-1].neggroup
            inposplusgroup = list(self.subclades[len(self.subclades)-1].neggroup)    
        for kit in self.plusgroup:
            inposplusgroup.append(kit)
        valid = False 
        if maxcount == 0:
            maxcount = len(inposplusgroup)-1
    
        #if maxcount == 0:
        #    self.foundAllSubclades = True
        #    self.failedQualityCheck = True
        #    return            
        (max_frequency_SNPs, maxcount) = max_frequency_SNP_list(ingroup, inposplusgroup, dfnotneg, maxcount, inconsistentSNPs, skipSNPs, allKits)
        #if maxcount == 0:
        #    self.foundAllSubclades = True
        #    return
         
        if len(max_frequency_SNPs) == 0:
            return -1
            
        newclade = SNPequivalents()
        add_SNPs_to_clade(newclade, max_frequency_SNPs, allSNPs, SNPnamesAll)
        print("{0} : {1}".format(newclade.name, maxcount) )

        (posplusgroup, neggroup, negplusgroup, posgroup, plusgroup) = group_kits(ingroup, newclade.primarySNP.position, df)
        newclade.addPosGroup(posgroup)            
        #if maxcount == -1: #singletons
        #    if len(negplusgroup) == 0:
        #        self.foundAllSubclades = True
        #        self.addSubclade(newclade)
        #        return
        newclade.addPlusGroup(plusgroup)
        newclade.addNegGroup(neggroup)
        return newclade    
           
    def checkIfAPosBNeg(self, SNPA, SNPB, df, dfpos):    
        hasSharedPos = False
    
        for kit in df.columns: 
            if df.ix[SNPA, kit] == 1 and df.ix[SNPB, kit] == 1:
                hasSharedPos = True
                break
        if hasSharedPos:        
            for kit in df.columns:          
                if SNPB == 13973248 and df.ix[SNPA, kit] == 1:
                    print(df.ix[SNPA, kit])
                    print(df.ix[SNPB, kit])
                if df.ix[SNPA, kit] == 1 and df.ix[SNPB, kit] == -1:
                    return True #bad SNP
        else:
            return False #Good SNP (brother clade)
        return False #Good SNP
           
           
             

    def findSubclade(self, df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll, dfpos):
        topskipSNPs = []
        valid = False
        while not valid:
            skipSNPs = list(topskipSNPs)
            newclade = self.findSubcladeCandidate( df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll, topskipSNPs, 0)
            if newclade == -1:
                return -1
            #Add the found cladecandidate to skipSNPs, and find the next candidate:        
            skipSNPs.append(newclade.primarySNP.position) 
            maxcount = len(newclade.posgroup) + len(newclade.plusgroup)      
            
            failedTest = False
            while not failedTest:
                        
                tempclade = self.findSubcladeCandidate(df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll, skipSNPs, maxcount)
                if tempclade == -1:
                    failedTest = True
                    break
                for SNP in tempclade.SNPs:
                    if tempclade.checkIfAPosBNeg(SNP.position, newclade.primarySNP.position, df, dfpos): #newclade failed
                        failedTest = True
                        break
                    skipSNPs.append(SNP.position)
                maxcount = len(tempclade.posgroup) + len(tempclade.plusgroup)
                if maxcount< len(newclade.posgroup) or maxcount == 1:
                    break
                 
            if not failedTest:
                while len(newclade.plusgroup)> 0:
                    
                    if newclade.findSubclade(df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll, dfpos) == -1: #no more subclades
                        for kit in newclade.plusgroup:
                            newclade.neggroup.append(kit)
                        break
                    newclade.fixPlusgroup()
                
                self.addSubclade(newclade) 
                valid = True
                print("{0} Added to {1}".format(newclade.name,self.name))
                            
            else:
                topskipSNPs.append(newclade.primarySNP.position)
                print("{0} failed".format(newclade.name))           
                
            
    def checkIfPositive(self,kit):
        if kit in self.posgroup:
            return 1
        else:
            sumofPos = 0
            for SNP in self.SNPs:
                if df.ix[SNP.position, kit] == 1:
                    return 1
            return 0
    def addToGroup(self, kit):
        notfinished = 1
    def addToPosGroup(self,kit):
        if kit not in self.posgroup:
            self.posgroup.append(kit)
            if len(self.subclades) > 0:
                for clade in self.subclades:
                    if clade.checkIfPositive == 1:
                        clade.addToPosGroup(kit)
                    else:
                        clade.addToGroup(kit)
    def movePlusToPosgroup(self,kit):
        if kit in self.plusgroup:
            self.posgroup.append(kit)
            self.plusgroup.remove(kit)
            if len(self.subclades) > 0:
                for clade in self.subclades:
                    if clade.checkIfPositive(kit) == 1:
                        clade.addToPosGroup(kit)
                    else:
                        clade.addToGroup(kit)                
    def fixPlusgroup(self):
        if len(self.plusgroup) > 0:
            plusgroupcopy = list(self.plusgroup)
            for kit in plusgroupcopy:
                for clade in self.subclades:
                    for SNP in clade.SNPs:
                        if df.ix[SNP.position, kit] == 1:
                            self.movePlusToPosgroup(kit)
        if len(self.plusgroup) > 0:
            return -1
    def checkifvalid(self, df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll):
        'checks if neggroup has positive in subclades'
        if len(self.plusgroup) > 0:
            while not self.foundAllSubclades:
                self.findSubclade(df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll)
            if self.failedQualityCheck:
                return False
            self.fixPlusgroup()
        return True        


allSNPs = []

with open('allSNPs.csv', 'r') as allSNPs_fp:

    SNPreader = csv.reader(allSNPs_fp, delimiter=',')
    for row in SNPreader:
        allSNPs.append(row)

#M = np.load('allSNPs.npy')
df = pd.read_pickle('allSNPs.pkl')

#print(df.ix[2650033, '333649'])
#print(df['333649'].head('8796078'))
#print("snp_map[1] = {0}".format(snp_map[1]))
#print("M[1,1] = {0}".format(M[1,1]))

def group_kits(startgroup, SNP, df):
    'splits the kit group into two groups depending on value of SNP'
    posgroup = list()
    posplusgroup = list()
    neggroup = list()
    negplusgroup = list()
    plusgroup = list()
    for kit in startgroup:
        if df.ix[SNP, kit] != -1:
            posplusgroup.append(kit)
        else:
            neggroup.append(kit)            
        if df.ix[SNP, kit] == 1:
            posgroup.append(kit)
        else:
            negplusgroup.append(kit) 
        if df.ix[SNP, kit] == 0 or np.isnan(df.ix[SNP, kit]):
            plusgroup.append(kit)
        
    
    return (posplusgroup, neggroup, negplusgroup, posgroup, plusgroup)

def frequency_SNP_list(kitgroup, df, maxcount = 0):
    sortedSNPlist = df[df == 1].loc[:,kitgroup].count(1).sort_values(ascending=False)
    if maxcount == 0:
        maxcount = sortedSNPlist.iloc[0]+1
    newsplitSNPs= sortedSNPlist[sortedSNPlist <= maxcount].index  
    for SNP in sortedSNPlist.index:
        if sortedSNPlist.loc[SNP] < maxcount:
            nextmaxcount = sortedSNPlist.loc[SNP]
            break
    return (newsplitSNPs, nextmaxcount)

def selectEquivalentSNPs(approved_SNPs, kitgroup, dfnotneg):
    mainSNP = approved_SNPs[0]
    out_SNPs = []
    for SNP in approved_SNPs:
        same = True
        for kit in kitgroup:
            if not dfnotneg.ix[SNP,kit] == dfnotneg.ix[mainSNP,kit]:
                same = False
                continue
        if same:
            out_SNPs.append(SNP)
        
    return out_SNPs
    
    
def max_frequency_SNP_list(kitgroup, posplusgroup, dfnotneg, maxcount, inconsistentSNPs, skipSNPs, allKits):
    max_frequency_SNPs = []    
    if len(posplusgroup) == 1:
        restgroup = list(allKits)
        #print(restgroup)
        restgroup.remove(posplusgroup[0])
        max_frequency_SNPs = find_exclusive_SNPs(posplusgroup, restgroup, df, dfpos, dfnotneg)
    else:
        sortedSNPlist = dfnotneg.loc[:,kitgroup].sum(1).sort_values(ascending=False)
        max_frequency_SNPs = sortedSNPlist[sortedSNPlist == maxcount].index 
        maxcount = min(sortedSNPlist.iloc[0], maxcount)
    if maxcount < 1:
        return ([],0)    
    
    approved_SNPs = []
    while len(approved_SNPs) == 0:
         

        for SNP in max_frequency_SNPs:
            if SNP == 9966805:
                pause = True            
                #print(dfpos.loc[SNP,kitgroup].count()>maxcount-20)
                #print(dfpos.loc[SNP,:].count())
                #print(dfpos.loc[SNP,posplusgroup].count())
            if SNP not in inconsistentSNPs and SNP not in skipSNPs and dfpos.loc[SNP,kitgroup].count()>maxcount-20 and dfpos.loc[SNP,:].count() == dfpos.loc[SNP,posplusgroup].count():
                approved_SNPs.append(SNP)
        if len(approved_SNPs) == 0:
            maxcount -= 1
        if maxcount < 1:
            return ([],0)
        max_frequency_SNPs = sortedSNPlist[sortedSNPlist == maxcount].index 
    out_SNPs = selectEquivalentSNPs(approved_SNPs, kitgroup, dfnotneg)
    #print(dfnotneg.loc[out_SNPs[0],kitgroup])
    if dfpos.loc[out_SNPs[0],:].count() == 1 and dfnotneg.loc[out_SNPs[0],kitgroup].count() == 1:
        maxcount = -1
    return (out_SNPs, maxcount)    


def find_first_consistent(newsplitSNPs, posgroup, dfpos, badSNPlist):
    success = False
    nextsplitSNP = False
    for SNP in newsplitSNPs:
        #if SNP == 14991735:
        #    h=0
        if SNP in badSNPlist:
            continue
        no_of_positives_in_all_kits = dfpos.loc[SNP,:].count()
        notposgroup  = list()   
        if no_of_positives_in_all_kits<300:
            for kit in dfpos.columns:
                if kit not in posgroup:
                    if dfpos.loc[SNP,kit] == 1:
                        y = 0
            
        
        no_of_positives_in_group = dfpos.loc[SNP,posgroup].count()    
        
        if no_of_positives_in_all_kits <= no_of_positives_in_group:
            success = True
            nextsplitSNP = SNP
            nextmaxcount = no_of_positives_in_group
            break
    return (nextsplitSNP, no_of_positives_in_group)

def find_exclusive_SNPs(posgroup, neggroup, df, dfpos, dfnotneg):
     
    return dfpos[dfpos.loc[:,neggroup].count(1) == 0].loc[:,posgroup].index
    

def add_SNPs_to_clade(clade, SNPlist, allSNPs, SNPnamesAll):
    for SNPposition in SNPlist:
        SNP = SNPPrototype(SNPposition)
        SNP.findNames(SNPnamesAll)
        allSNPs.append(SNP)
        clade.addSNP(SNP)
        if SNP.position in preferredNames:
            SNP.setPrimaryName(preferredNames[SNP.position])
            clade.setPrimarySNP(SNP) 
    


def get_name(SNP, SNPnames):
    if str(SNP) in SNPnames:
        return SNPnames[str(SNP)]
    else:
        return SNP

def find_largest_child_SNP(posplusgroup, df, dfpos, maxcount, badSNPlist):
    neggroup = []
    maxcountin = maxcount
    while len(neggroup) == 0:
        (newsplitSNPs, maxcount) = frequency_SNP_list(posplusgroup, df, maxcount)
        if maxcount == 1:
            return ([],[],[],[],False,maxcountin,[])
        (newsplitSNP, blamaxcount) = find_first_consistent(newsplitSNPs, posplusgroup, dfpos, badSNPlist)
        if blamaxcount == 1:
            return ([],[],[],[],False,maxcountin,[])
        (posplusgroup, neggroup, negplusgroup, posgroup, plusgroup) = group_kits(posplusgroup, newsplitSNP, df)
    return (posplusgroup, neggroup, negplusgroup, posgroup, newsplitSNP, maxcount, plusgroup)
    
def find_next_child(ingroup, df, dfpos, maxcount, badSNPlist):
    neggroup = []
    if len(ingroup)>1:
        while len(neggroup) == 0:
            (posplusgroup, neggroup, negplusgroup, posgroup, newsplitSNP, maxcount, plusgroup) =find_largest_child_SNP(ingroup, df, dfpos, maxcount, badSNPlist)
            if maxcount == 0:
                return (ingroup, neggroup, neggroup, False, 0)
    else:
        return (ingroup, neggroup, neggroup, False, 0)
    return (posplusgroup, neggroup, negplusgroup, newsplitSNP, maxcount)

def get_all_children_old(posplusgroup, df, dfpos, nextmaxcount, badSNPlist, parentname, SNPnames):
    #largest child:
    if len(posplusgroup) > 1:
        neggroup = []
        while len(neggroup) == 0:
            (posplusgroup, neggroup, negplusgroup, posgroup, newsplitSNP, nextmaxcount, plusgroup) = find_largest_child_SNP(posplusgroup, df, dfpos, nextmaxcount, badSNPlist)
            if nextmaxcount == 0:
                return ([], [], [])
        newname = "{0}>{1}".format(parentname, get_name(newsplitSNP,SNPnames))
        lastSNP = newsplitSNP
        #print(newname)
        SNPchildren = [newsplitSNP]
        nameChildren = [newname]
        posplusgroupchildren = [posplusgroup]
        
        while len(negplusgroup) > 0:
            (posplusgroup, neggroup, negplusgroup, newsplitSNP, nextmaxcount) = find_next_child(negplusgroup, df, dfpos, nextmaxcount, badSNPlist)        
            if lastSNP != newsplitSNP and len(negplusgroup) > 0 and newsplitSNP != False:
                #print(neggroup)
                newname = "{0}>{1}".format(parentname, get_name(newsplitSNP,SNPnames))
                lastSNP = newsplitSNP
                #print(newname)    
                SNPchildren.append(newsplitSNP)
                nameChildren.append(newname)
                posplusgroupchildren.append(posplusgroup)
        return (SNPchildren, nameChildren, posplusgroupchildren)
    else: return ([], [], [])

def get_all_children(ingroup, df, dfpos, maxcount):
    if len(ingroup)>1:
        neggroup = ["dummy"]
        childSNPs = list()
        while len(neggroup)>0:
            (posplusgroup, neggroup, negplusgroup, posgroup, newsplitSNP, maxcount, plusgroup) = find_largest_child_SNP(ingroup, df, dfpos, maxcount, badSNPlist)
            if newsplitSNP == False: #only single children
                (newsplitSNPs, maxcount2) = frequency_SNP_list(ingroup, df, maxcount)
                (newsplitSNP, blamaxcount) = find_first_consistent(newsplitSNPs, ingroup, dfpos, badSNPlist)  
                if blamaxcount > 1:
                    #print(get_name(newsplitSNP, SNPnames))
                    childSNPs.append(get_name(newsplitSNP, SNPnames))
                break;
            #print(get_name(newsplitSNP, SNPnames))
            childSNPs.append(get_name(newsplitSNP, SNPnames))
            if len(plusgroup) > 0: #is it nessecary to check where plusgroup belongs
                
                not_in_children = is_plusgroup_in_children(posplusgroup, df, dfpos, maxcount, badSNPlist, plusgroup, SNPnames, df)
                for kit in plusgroup:
                    if kit in not_in_children:
                        neggroup.append(kit)
                    else:
                        posgroup.append(kit)
            ingroup = neggroup
        return (childSNPs)
    else: return ([])
                    

def is_plusgroup_in_children(ingroup, dfred, dfredpos, maxcount, badSNPlist, inplusgroup, SNPnames, df): 
    #remove inplusgroup from df:
    dfreduced = dfred 
    for pluskit in inplusgroup:
        
        if pluskit in dfred.columns:
            dfreduced = dfreduced.drop(pluskit,1)  
          
    dfreducedpos = dfreduced[dfreduced == 1]  
    out = list()
    for pluskit in inplusgroup:
        print(pluskit)
        analysisgroup = list()
        for kit in ingroup:
            if kit not in inplusgroup:
                analysisgroup.append(kit)
        if len(analysisgroup)<2: #Actually do direct comparison.
            if pluskit not in out:
                out.append(pluskit)
            continue
        (posplusgroup, neggroup, negplusgroup, posgroup, newsplitSNP, nextmaxcount, plusgroup) = find_largest_child_SNP(analysisgroup, dfreduced, dfreducedpos, maxcount, badSNPlist)
        
        if newsplitSNP == False:
            if pluskit not in out:
                out.append(pluskit)
            continue
        a = df.ix[newsplitSNP, pluskit]
        x = get_name(newsplitSNP, SNPnames)
        print(x)
        if df.ix[newsplitSNP, pluskit] != 1:  #option A - remove plusses V
            
            nowinneggroup = False
            nowinplusgroup = False
            if df.ix[newsplitSNP, pluskit] == -1:
                nowinneggroup= True   #option b - test neggroup  V           
            if df.ix[newsplitSNP, pluskit] == 0 or np.isnan(df.ix[newsplitSNP, pluskit]):
                nowinplusgroup = True   #option c - test children/neg        
            
            #find next child with neggroup2 - check status of nowinneggroup2
            if nowinneggroup:
                print("neg")
                for kit in is_plusgroup_in_children(neggroup, dfreduced, dfreducedpos, nextmaxcount, badSNPlist, [pluskit], SNPnames, df):
                    if kit not in out:
                        out.append(kit)
            if nowinplusgroup:
                print("plus")
                stillinplusgroup = False
                for kit in is_plusgroup_in_children(neggroup, dfreduced, dfreducedpos, nextmaxcount, badSNPlist, [pluskit], SNPnames, df):
                
                    for kit in is_plusgroup_in_children(posplusgroup, dfreduced, dfreducedpos, nextmaxcount, badSNPlist, [pluskit], SNPnames, df):
                        if kit not in out:
                            out.append(kit)                
    return out


SNPnames = dict()
SNPnamesAll = dict()

with open('SNPNames.csv', 'r') as SNPnames_fp:

    Namesreader = csv.reader(SNPnames_fp, delimiter=',')
    for row in Namesreader:
        SNPnames[row[1]] = row[0]
        if row[1] in SNPnamesAll:
            SNPnamesAll[row[1]].append(row[0])
        else:
            SNPnamesAll[row[1]] = [row[0]]

#print(SNPnamesAll)
#no_of_positives = defaultdict(int)

#for SNP in allSNPs:
    #if SNP[2] == "pos":
        #no_of_positives[SNP[1]] += 1

#sorted_no_of_positives = sorted(no_of_positives.items(), key=operator.itemgetter(1), reverse=True)


#out = []
#for SNP, count in sorted_no_of_positives:
    #if str(SNP) in SNPnames:
        #out.append([SNPnames[str(SNP)],count])
    #else:
        #out.append([SNP, count])

#with open('frequencysortednames.csv', 'w') as fp:
    #writer = csv.writer(fp)
    #writer.writerows(out)

#maxpos = max(no_of_positives.values())
#print(maxpos)

inconsistentSNPs = [22487098, 13211699, 13220892, 13238659, 24447989, 22269795, 13204326, 18060384, 18060417, 15172515, 9031181, 14577177, 14993358, 14399057, 3714315, 18060390, 14399074, 14399077, 22346168, 13239873, 14367269, 13204380, 28788695, 22319170, 22319171, 14399061, 15080010, 14394512, 14956117, 14753989, 24453842, 14989721, 22318978, 14639427, 7321330, 22246199, 22318641, 22318642, 22318643, 23354909, 22269796, 22444659, 14398980, 14399063, 15094399, 15094400, 22293981, 28787127, 28787141, 22233732, 22262889, 14415575, 22246045, 14825721, 14071351, 13687378, 13204363, 9227801, 22259729, 28788643, 22235226, 9316833, 15027529, 22229722, 14624254, 15202608, 9810794, 14424045, 22269664, 14392996, 9992071, 22319065, 22319066, 9296025, 14646338, 9374685, 23813856, 14478992, 14479010, 14478968, 14478959, 14479024, 14610658, 14479048, 14071348, 26123300, 24381403, 22436085, 14479072, 13707988, 17860729, 9239350, 9852981, 14926420, 14630342, 14630334, 14630357, 22240888, 14479079, 14624294, 14423856, 14981750, 16607930, 22316451, 22436108, 23311208, 13271933, 22261622, 13494176, 7717685, 9981785, 14561760, 14598808, 14479100, 23673867, 28802526, 14933671, 9350171, 9937995, 15558499, 23897432, 17364720, 28792932, 23488865, 17467411, 17467482, 17467467, 19963306, 9311208, 13521711, 9311383, 13204205, 13717500, 22288897, 28804383, 28818135, 14435782, 14435790, 14435792, 14435818, 13268521, 22419921, 13142458, 7307494, 13142438, 13142443, 13686963, 13142433, 28587749, 7534406, 28793052, 13204196, 28587738, 13687415, 24475333, 22301589, 7717741, 14480546, 19048394, 13502130, 22304778, 19936492, 22316647, 28816999, 28817000, 17275007, 7164283, 7730475, 22241361, 22302276, 22302653, 19436738, 7643546, 10012012, 13851646, 13687433, 13687807, 13687826, 6933361, 13312756, 9321336, 9981898, 16692107, 22299085, 22318301, 22299728, 22303133, 13687570, 13851655, 22232987, 22233360, 22299844, 22299739, 21153096, 10012071, 13804879, 13687977, 28816948, 16692095, 16494396, 22234486, 22230508, 22235382, 22310538, 13294106, 13804809, 16494398, 22256872, 26124474, 22420062, 22291623, 22435876, 13294119, 10012098, 16692088, 22299750, 19376214, 13851668, 15799233, 19968557, 22299752, 18892281, 23709848, 6141754, 10012109, 22240008, 14242418, 16505988, 14242413, 10012140, 10012141, 28793049, 13804902, 19048392, 22232260, 25931596, 15558495, 2824726, 9231864, 16818679, 19469877, 22490457, 22420067, 5690586, 13687497, 22239394, 18196247, 28817559, 13268577, 16692081, 13573108, 10012138, 19947053, 22487462, 18892280, 7157050, 22229844, 22319609, 17054892, 9320413, 13301159, 13301165, 22259014, 23311204, 13866490, 13804788, 27782617, 14436052, 16692103, 22299777, 15278826, 22487541, 22229834, 22309944, 19103890, 24961157, 14242408, 23768775, 22435870, 14102213, 22316424, 19377439, 6933386, 13722195, 13323548, 28802511, 18687682, 14926489, 22299774, 7307496, 18149594, 14298753, 9220362, 15226343, 16441971, 22420085, 17946919, 24278193, 13573152, 13435187, 8479424, 9207152, 9320093, 7136939, 22288750, 9994759, 13301185, 17274963, 15320664, 25080669, 13664772, 18149590, 22469341, 13294025, 13307425, 20000369, 26134809, 19841099, 22229828, 19971812, 22460839, 21477232, 6753258, 21152868, 16673063, 14612163, 18982587, 18982595, 22469846, 22425646, 22425924, 23136653, 13317375, 13833231, 13801582, 22281790, 22288793, 13728534, 59021925, 22438672, 13448404, 13680834, 21714748, 22802901, 22473179, 22473273, 22473287, 17386058, 13689668, 13851574, 9330881, 9372135, 13268593, 22477710, 22477716, 22477751, 22477758, 13302072, 22461574, 19978363, 18588420, 17794253, 15226371, 22478089, 22253868, 22253976, 22256943, 22336565, 15810962, 22467097, 22467126, 22467191, 22462613, 22462616, 22462651, 18589702, 15226328, 15226331, 15226335, 15226338, 19209917, 19209925, 19209943, 2888657, 22487366, 21152862, 13721754, 18687686, 26179808, 22353065, 13833192, 13833214, 28818131, 22281789, 13833226, 13833229, 13865916, 22436396, 18995237, 14706275, 16394317, 13717184, 13717866, 13711202, 13713916, 7020430, 13866440, 22439297, 22441176, 22441305, 15278816, 13664768, 13470087, 15491978, 9843195, 23832815, 23832820, 21138770, 22802902, 22473255, 9231987, 13481556, 13664774, 13664781, 13664819, 28755129, 9754793, 13142444, 28804016, 21285190, 22473949, 22474043, 13312446, 13312452, 13312469, 22272627, 13309262, 22460570, 22460598, 19986086, 21846452, 17192869, 22287932, 22287951, 16378861, 22259728, 19436768, 18518976, 18372055, 26282891, 21477226, 16494359, 13494230, 22318649, 14287205, 6105265, 18283909, 23311212, 22242160, 22462538, 22462549, 22464022, 15818409, 9730834, 8432962, 8568432, 15255773, 19706640, 15226705, 22287415, 22441018, 19423834, 13704254, 22251944, 14294241, 21262424, 8650100, 8419806, 21398583, 13869393, 17387296, 16441973, 9654393, 14609823, 17638126, 7717737, 17850011, 28802514, 21152967, 21152971, 17467526, 18119306, 27232111, 9745685, 22459353, 22459365, 22459368, 24071473, 8273411, 7136945, 13309174, 25303064, 22459019, 22459048, 22459063, 22459064, 22459098, 22459121, 22459127, 22460668, 22460746, 23037438, 22461050, 24553622, 22270062, 22270127, 22270162, 22436961, 17649187, 18419695, 15226253, 22458391, 15477748, 3544787, 24553539, 24600517, 13654538, 17065062, 22289157, 14491211, 22439953, 22450575, 13417250, 22254661, 22439997, 22473346, 22474638, 13410742, 23415312, 9179212, 16410277, 9442049]


preferredNames = {8796078 : "U106", 7246726 : "Z381", 15780341 : "Z156", 14991735 : "Z18"}

allSNPs = []
allKits = df.columns
dfnotneg = df != -1
#print(dfnotneg)
dfpos = df[df == 1]
(max_frequency_SNPs, maxcount) = max_frequency_SNP_list(allKits, allKits, dfnotneg, 100000, inconsistentSNPs, [], allKits)

#print(len(allKits))
#print(max_frequency_SNPs)
topclade = SNPequivalents() #All the SNPs that are shared by everyone
add_SNPs_to_clade(topclade, max_frequency_SNPs, allSNPs, SNPnamesAll)

topclade.addPosGroup(allKits)

#find all the ones with one negative. Most are false negatives, except for those that are onle negative for kit 87721
(max_frequency_SNPs, maxcount) = max_frequency_SNP_list(allKits, allKits, dfnotneg, maxcount-1, inconsistentSNPs, [], allKits) 


#for SNP in max_frequency_SNPs:
#    (posplusgroup, neggroup, negplusgroup, posgroup, plusgroup) = group_kits(allKits, SNP, df)
#    print("{0}:{1},{2}".format(SNP,neggroup,plusgroup))
   
newclade = SNPequivalents()
#We ignore all other SNPs than 14641193 and 8796078
U106equiv = [14641193,8796078]
add_SNPs_to_clade(newclade, U106equiv, allSNPs, SNPnamesAll)

        
(posplusgroup, neggroup, negplusgroup, posgroup, plusgroup) = group_kits(allKits, newclade.primarySNP.position, df)
newclade.addPosGroup(posgroup)
newclade.addPlusGroup(plusgroup)
newclade.addNegGroup(neggroup)

topclade.addSubclade(newclade)

#get singeltons of the U106 kit
singletons = find_exclusive_SNPs(neggroup, posplusgroup, df, dfpos, dfnotneg)
newclade = SNPequivalents()
add_SNPs_to_clade(newclade, singletons, allSNPs, SNPnamesAll)
newclade.addPosGroup(neggroup)
topclade.addSubclade(newclade)

# now we have a topclade with U106 and one lone kit as subclades.

# investigate U106:
clade = topclade.subclades[0]
group = clade.posgroup

#Inconsistent SNPs form csv-sheet. should be autogenerated somehow in the future
clade.findSubclade( df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll, dfpos)

while len(clade.subclades[len(clade.subclades)-1].neggroup) >0:
    clade.findSubclade( df, dfnotneg, inconsistentSNPs, allSNPs, SNPnamesAll, dfpos)

               


badSNPlist = [14624254,18196247]

startSNP = 15780341 #U106 = 8796078, Z381 = 7246726
startName = "Z156"

(posplusgroup, neggroup, negplusgroup, posgroup, plusgroup) = group_kits(df.columns, startSNP, df) #get groups from startSNP - neg is rest of kits
(newsplitSNPs, nextmaxcount1) = frequency_SNP_list(posplusgroup, df, 0)
(newsplitSNP, blamaxcount) = find_first_consistent(newsplitSNPs, posplusgroup, dfpos, badSNPlist) #returns startSNP (or equivalent)
(posplusgroup1, neggroup1, negplusgroup1, posgroup1, plusgroup1) = group_kits(posplusgroup, newsplitSNP, df) #get groups from startSNP - neg is empty
ingroup = posplusgroup1
(childSNPs) = get_all_children(ingroup, df, dfpos, nextmaxcount1)
print(childSNPs)

sys.exit("hi")
neggroup1 = ["dummy"]
while len(neggroup1)>0:
    (posplusgroup1, neggroup1, negplusgroup1, posgroup1, newsplitSNP1, nextmaxcount1, plusgroup1) = find_largest_child_SNP(ingroup, df, dfpos, nextmaxcount1, badSNPlist)
    if newsplitSNP1 == False: #only single children
        (newsplitSNPs, maxcount) = frequency_SNP_list(ingroup, df, nextmaxcount1)
        (newsplitSNP, blamaxcount) = find_first_consistent(newsplitSNPs, ingroup, dfpos, badSNPlist)  
        if blamaxcount > 1:
            print(get_name(newsplitSNP, SNPnames))
        break;
    print(get_name(newsplitSNP1, SNPnames))
    if len(plusgroup1) > 0: #is it nessecary to check where plusgroup belongs
        
        not_in_children = is_plusgroup_in_children(posplusgroup1, df, dfpos, nextmaxcount1, badSNPlist, plusgroup1, SNPnames, df)
        for kit in plusgroup1:
            if kit in not_in_children:
                neggroup1.append(kit)
            else:
                posgroup1.append(kit)
    ingroup = neggroup1            
       
  
  
    
sys.exit("Stop") 
if a ==0:
    (posplusgroup2, neggroup2, negplusgroup2, posgroup2, newsplitSNP2, nextmaxcount2, plusgroup2) = find_largest_child_SNP(posplusgroup1, df, dfpos, nextmaxcount1, badSNPlist)
    stillinplusgroup1 = list()
    for kit in plusgroup1:
        if kit not in posgroup2:  #option A - remove plusses
            stillinplusgroup1.append(kit)
    if len(stillinplusgroup1) > 0:
        
        nowinneggroup2 = list()
        for kit in stillinplusgroup1:
            if kit in neggroup2:
                nowinneggroup2.append(kit)   #option b - test neggroup             
            if kit in plusgroup2:
                nowinplusgroup2.append(kit)   #option c - test children/neg        
        #find next child with neggroup2 - check status of nowinneggroup2
        if len(nowinneggroup2) >0: 
            (posplusgroup3, neggroup3, negplusgroup3, posgroup3, newsplitSNP3, nextmaxcount3, plusgroup3) = find_largest_child_SNP(neggroup2, df, dfpos, nextmaxcount2, badSNPlist)
            stillinplusgroup3 = list()
            for kit in nowinneggroup2:
                if kit not in posgroup3:
                    stillinplusgroup3.append(kit) #option A remove plusses
                    if len(stillinplusgroup3) > 0:
                        nowinneggroup3 = list()
                        nowinplusgroup3  = LIST()
                        for kit in stillinplusgroup3:
                            if kit in neggroup3:
                                nowinneggroup3.append(kit)   #option b - test neggroup             
                            if kit in plusgroup3:
                                nowinplusgroup3.append(kit)   #option c - test children/neg        
                        #find next child with neggroup2 - check status of nowinneggroup2        
                
                
#if len(nowinplusgroup2) == 0:
#    success = True
#else

        
else:
#find next child with neggroup1
    h = 0
            
            



startSNP = 13221656 #U106 = 8796078, Z381 = 7246726
startName = "Z160"
dfpos = df[df == 1]
#Get U106:
(posplusgroup, neggroup, negplusgroup, posgroup, plusgroup) = group_kits(df.columns, startSNP, df)
print(get_name(startSNP,SNPnames))

(newsplitSNPs, nextmaxcount) = frequency_SNP_list(posplusgroup, df, 0)
(newsplitSNP, blamaxcount) = find_first_consistent(newsplitSNPs, posplusgroup, dfpos, badSNPlist)
(posplusgroup, neggroup, negplusgroup, posgroup, plusgroup) = group_kits(posplusgroup, newsplitSNP, df)

fullname = get_name(startSNP,SNPnames)
lastSNP = get_name(startSNP,SNPnames)
#levels = (SNPchildren, nameChildren, posplusgroupchildren)
levels = []
levels.append([[startSNP], [startName], [posplusgroup]])
i = 1
newSplits = 1
while newSplits > 0:
    newSplits = 0
    levels.append([])
    levels[i].append([])
    levels[i].append([])
    levels[i].append([])
            
    for j in range(len(levels[i-1][0])):    
        (SNPchildren, nameChildren, posplusgroupchildren) = get_all_children(levels[i-1][2][j], df, dfpos, 0, badSNPlist, levels[i-1][1][j], SNPnames)
        if nameChildren == ['U106>18645737>9761816']:
            pause = True
        if len(SNPchildren) > 0:
            for SNP in SNPchildren:            
                levels[i][0].append( SNP)
            for name in nameChildren:   
                levels[i][1].append( name)
            for posplusgroup in posplusgroupchildren:
                levels[i][2].append( posplusgroup)
            print(nameChildren)
            newSplits += 1
    i = i+1
        
print("done")
#level2 = get_all_children(posplusgroup, df, dfpos, 0, badSNPlist, fullname, SNPnames)
#print(level2[1])
#level3 =[]
#for i in range(len(level2[0])):
    #level3.append(get_all_children(level2[2][i], df, dfpos, 0, badSNPlist, level2[1][i], SNPnames))
    #print(level3[i][1])
    
    
    
    
    
    
#while len(neggroup) == 0:
    #(posplusgroup, neggroup, negplusgroup, newsplitSNP, nextmaxcount) = find_largest_child_SNP(posplusgroup, df, dfpos, nextmaxcount, badSNPlist)
#if lastSNP != newsplitSNP and len(neggroup) > 0:
        #print(neggroup)
        #newname = "{0}>{1}".format(fullname, get_name(newsplitSNP,SNPnames))
        #lastSNP = newsplitSNP
        #print(newname)


#while len(negplusgroup) > 0:
    #(posplusgroup, neggroup, negplusgroup, newsplitSNP, nextmaxcount) = find_next_child(negplusgroup, df, dfpos, nextmaxcount, badSNPlist)        
    #if lastSNP != newsplitSNP and len(neggroup) > 0:
        #print(neggroup)
        #newname = "{0}>{1}".format(fullname, get_name(newsplitSNP,SNPnames))
        #lastSNP = newsplitSNP
        #print(newname)        
    
    
    
#while nextmaxcount >0 :
#    (posplusgroup, neggroup, newsplitSNP, nextmaxcount) = find_largest_child_SNP(posplusgroup, df, dfpos, nextmaxcount, badSNPlist)
#    if lastSNP != newsplitSNP and len(neggroup) > 0:
        #print(neggroup)
#        fullname = "{0}>{1}".format(fullname, get_name(newsplitSNP,SNPnames))
#        lastSNP = newsplitSNP
#        print(fullname)
        
      
    

#if newsplitSNP == False:
 #   (newsplitSNPs, maxcount) = frequency_SNP_list(posgroup, df, maxcount)
  #  print(newsplitSNPs)
   # newsplitSNP = find_first_consistent(newsplitSNPs, posgroup, df)
    #print(newsplitSNP)

#print(newsplitSNPs)




#(posgroup, neggroup) = group_kits(posgroup, newsplitSNP, df)
#print(posgroup)