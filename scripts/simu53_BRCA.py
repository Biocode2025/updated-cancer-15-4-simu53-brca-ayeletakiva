#functions:
import random

#פונקציה המקבלת קובץ דנ"א והופכת אותו לרצף רנ"א
def DNA_RNA_Cod(DNA):
  delimiter = ""
  line_list = []
  DNA = DNA.upper()
  for ch in DNA:
    if ch == "T":
      line_list.append("U")
    else:
      line_list.append(ch)
  RNA = delimiter.join(line_list)
  return RNA

#פונקציה שמקבלת קובץ חומצות אמיניות ואת הקיצורים שלהן ומכניסה אותן למילון גלובלי
def Read_dict():
  global RNA_codon_table
  file = open("data/codon_AA.txt")
  for line in file:
    line = line.rstrip("\r\n")
    (codon,dev,AA) = line.partition("\t")
    RNA_codon_table[codon] = AA
  file.close()

#פונקציה המקבלת רצף רנ"א ומתרגמת אותו לחומצות האמיניות מהן הוא מורכב
def RNA_prot(RNA):
  AA_protein = ""
  for i in range(0, len(RNA), 3):
    codon = RNA[i:i+3]
    if len(codon) == 3:
      if codon in RNA_codon_table:
        if RNA_codon_table[codon] == "*":
          AA_protein += (RNA_codon_table[codon])
          break
        else:
          AA_protein += (RNA_codon_table[codon])
  return AA_protein

#פונקציה המקבלת את רצף הדנ"א, בוחרת מיקום אקראי ברצף ומחליפה בו את הנוקלאוטיד באופן רנדומלי, מחזירה את רצף הגנום עם המוטציה
def Mutate_DNA(seq):
  base_list = ["A", "T", "C", "G"]
  rand_base = random.randrange(0,len(seq))
  if seq[rand_base] == "A":
    base_list.remove("A")
  elif seq[rand_base] == "T":
    base_list.remove("T")
  elif seq[rand_base] == "C":
    base_list.remove("C")
  elif seq[rand_base] == "G":
    base_list.remove("G")
  new_base = random.choice(base_list)
  mut_seq = seq[0:rand_base] + new_base + seq[rand_base+1:]
  return mut_seq

# פונקציה המסירה נוקלאוטיד אחד עד שלושה במקום אקראי ברצף
def Delete_DNA(seq):
  nuc_amount = random.randrange(1,4)
  rand_base = random.randrange(0,len(seq))
  mut_seq = seq[0:rand_base] + seq[rand_base+nuc_amount:]
  return mut_seq

#פונקציה המכניסה במיקום אקראי לרצף נוקלאוטיד אחד עד שלושה נוספים
def Insert_DNA(seq):
  base_list = ["A", "T", "C", "G"]
  rand_base = random.randrange(0,len(seq))
  nuc_amount = random.randrange(1,4)
  new_bases = ""
  for i in range(nuc_amount):
    new_base = random.choice(base_list)
    new_bases += new_base
  mut_seq = seq[0:rand_base] + new_bases + seq[rand_base:]
  return mut_seq

#פונקציה המשווה ובודקת כמה הבדלים יש בין שני רצפים ומחזירה את מספר ההבדלים הקיימים ביניהם
def Comp_seq(old,new):
  diff = 0
  min_length = min(len(old), len(new))
  for i in range(min_length):
    if old[i] != new[i]:
      diff += 1
    else:
      continue
  diff += abs(len(old) - len(new))
  return diff

### main program ###

#יצירת מילון הקודונים
RNA_codon_table = {}
Read_dict()

#הכנסת הקובץ לשורת סטרינג אחת
file = open("data/human_p53_coding.txt")
p53 = ""
for line in file:
  if line[0] == ">":
    continue
  else:
    line = line.upper()
    line = line.rstrip("\r\n")
    p53 += line

times = 1000
gene_occ = input("Does the Female has a BRCA1,2 mutation? (Y=Yes, N=No)")

#הוספת מוטציה רנדומלית
og_RNA = DNA_RNA_Cod(p53)
og_prot = RNA_prot(og_RNA)
tot_day_count = 0
#במקרה ויש מוטציה אחת
if gene_occ == "yes":
  for i in range(times):
    m_p53 = p53
    day_count = 0
    diff = 0
    while diff == 0:
      day_count += 1
      chance_mut = random.randint(1,10000)
      if chance_mut == 10000:
        chance = random.randint(1,100)
        if chance <= 98:
          m_p53 = Mutate_DNA(m_p53)
        elif chance == 99:
          m_p53 = Delete_DNA(m_p53)
        elif chance == 100:
          m_p53 = Insert_DNA(m_p53)
        mut_RNA = DNA_RNA_Cod(m_p53)
        mut_prot = RNA_prot(mut_RNA)
        diff = Comp_seq(og_prot, mut_prot)
        if diff > 0:
          break
      else:
        continue
    tot_day_count += day_count
  av_days = tot_day_count/times
  years = av_days/365
  print("For a female that does have BRCA1,2 Mutation: \nThe mutation that will change the P53 protein will take in average %d years." %(years))
  
#במקרה ואין בכלל מוטציה
elif gene_occ == "no":
  for i in range(times):
    m_p53 = p53
    day_count = 0
    diff = 0
    while diff <= 1:
      day_count += 1
      chance_mut = random.randint(1,10000)
      if chance_mut == 10000:
        chance = random.randint(1,100)
        if chance <= 98:
          m_p53 = Mutate_DNA(m_p53)
        elif chance == 99:
          m_p53 = Delete_DNA(m_p53)
        elif chance == 100:
          m_p53 = Insert_DNA(m_p53)
        mut_RNA = DNA_RNA_Cod(m_p53)
        mut_prot = RNA_prot(mut_RNA)
        diff = Comp_seq(og_prot, mut_prot)
        if diff == 2:
          break
      else:
        continue
    tot_day_count += day_count
  av_days = tot_day_count/times
  years = av_days/365
  print("For a female that does have BRCA1,2 Mutation: \nThe mutation that will change the P53 protein will take in average %d years." %(years))
