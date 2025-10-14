import csv
import copy
import sys
import os
import numpy as np
import copy
import random
import time
import math
import pickle
import shutil

random.seed(1)
np.random.seed(1)


def DA(intermed_dir, r_monte):
  '''
  New implementation of DA student optimal with strict preferences and scores
  '''

  # students.csv structure:
  # student_id
  # s1
  # s2
  # ...
  r_monte = int(r_monte)
  student_id_path = os.path.normpath(os.path.join(intermed_dir, f'students_py_{r_monte}.csv'))
  with open(student_id_path, 'r') as f:
    reader = csv.reader(f)
    next(reader)
    students = [row[0] for row in reader]

  #shutil.rmtree(student_id_path)
  

  # capacity.csv structure:
  # school_id, capacity
  # c1, 2
  # c2, 3
  # ...

  cap_path = os.path.normpath(os.path.join(intermed_dir, f'capacity_{r_monte}.csv'))
  cap = {}
  with open(cap_path, 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
      cap[row[0]] = int(row[1])

  #shutil.rmtree(cap_path)
  


  # rol.csv structure:
  # id, rank, option_id
  # s1, 1, c2
  # s1, 2, c1
  # c1, 1, s2
  # ...
  pref = {}
  rol_path = os.path.normpath(os.path.join(intermed_dir, f'rol_py_{r_monte}.csv'))
  with open(rol_path, 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
      if row[0] not in pref:
        pref[row[0]] = {}
      pref[row[0]][int(row[1])] = row[2]

  #shutil.rmtree(rol_path)
  
  # scores.csv structure:
  # student_id, score
  # s1, 85
  # s2, 88
  # ...

  scores = {}
  scores_path = os.path.normpath(os.path.join(intermed_dir, f'scores_py_{r_monte}.csv'))
  with open(scores_path, 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
      scores[row[0]] = float(row[1])
  # shutil.rmtree(scores_path)
  

 # Run the DA algorithm

 # Current proposal of students
  cp = {i: 1 for i in students}  

  # Match status dictionary
  match, is_matched = {i: None for i in students}, {i: False for i in students}

  # Each student proposes to its current top school
  proposals = {c: [] for c in cap} # Dictionary of school id and list of students who proposed to the school
  rejected = copy.copy(students)
  it = 0
  while True:
    it += 1
    print("")
    print(f"Proposal Round {it}: {len(rejected)} Unmatched")
    #print("Proposals that schools currently hold:")
    #print("Before new proposals.")
    # print([(c, len(proposals[c])) for c in cap])

    for m in rejected:
      if is_matched[m]:
        continue
      try:
        proposals[pref[m][cp[m]]].append(m)
      except:
        print(m, cp[m], pref[m][cp[m]])
        sys.exit(1)

    #print(f"After new proposals.")
    #print([(c, len(proposals[c])) for c in cap])

    rejected = []
    # Each school rejects those students who are not in the top based on scores
    for c in cap:
      # Sort students based on scores
      students_with_scores = [(student, scores[student])
                              for student in proposals[c]]  
      students_with_scores.sort(key=lambda x: x[1], reverse=True)
      top = [student for student, _ in students_with_scores]

      # if the number of proposals received is more than the vacancies
      if len(proposals[c]) > cap[c]:
        for i in range(cap[c], len(proposals[c])):
          proposals[c].remove(top[i])
          rejected.append(top[i])
    
    # The rejected students propose to their next top choice
    for m in rejected:
      cp[m] += 1
      if cp[m] not in pref[m]:
        is_matched[m] = True # Student is matched to the ouside option (None) 

    stop = all(is_matched[i] for i in students) # stop = True if all students are matched to the outside option
    if stop or len(rejected) == 0:
      break

  for m in students:
    if is_matched[m]:  # Student is matched to the ouside option (None), next iteration
      continue
    match[m], is_matched[m] = pref[m][cp[m]], True

  # Save the match to a CSV
  with open(os.path.normpath(os.path.join(intermed_dir, f'match_{r_monte}.csv')), 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['student_id', 'matched_school_id', 'matched_sch_type'])
    for student, school in match.items():
      school_id_actual = int(school.split("_")[0][1:])  
      writer.writerow([student, school_id_actual, school])

 

    
#if __name__ == "__main__":
#  intermed_dir = "Z:/decentralized_choice/data/intermediate/counterfactual_matching"
#  DA(intermed_dir, 1 )
