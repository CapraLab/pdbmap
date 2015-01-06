#!/usr/bin/python2.7
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=10000mb
#PBS -l walltime=5:00:00:00
import sys,os,csv,MySQLdb,math
os.chdir('/projects/Bush_eQTL/sivleyrm/projects/pdbmap/scripts')

# Adds GD to SNN and SD to GNN to the provided STEVE results file
# Uses gwar-dev pdbmap_v9
steve_finname = '../results/pdbmap-v9_steve_20140616-13/nearest_neighbors.txt'
steve_foutname = "%s.extended.txt"%'.'.join(steve_finname.split('.')[:-1])
con = MySQLdb.connect(host='gwar-dev.mc.vanderbilt.edu',user='mike',passwd='cheezburger',db='pdbmap_v9')
c = con.cursor()
fout = open(steve_foutname,'wb')
writer = csv.writer(fout,delimiter='\t')
with open(steve_finname,'rb') as fin:
  header = "%s\tGNN_structid\tGNN_biounit\tGNN_model\tGNN_x\tGNN_y\tGNN_z\tSD2GNN\tSNN_chr\tSNN_start\tGD2SNN\n"%fin.readline().strip()
  fout.write(header)
  reader = csv.reader(fin,delimiter='\t')
  for row in reader:
    VAR,GNN,SNN = row[0],row[2],row[7]
    if not SNN:
      gd2snn = 'NA'
    else:
      # Calculate GD to SNN
      c = con.cursor()
      query  = "SELECT chr,start FROM GenomicData "
      query += "WHERE label='1kg' AND name=%s LIMIT 1"
      c.execute(query,(VAR,))
      var_gloc = c.fetchone()
      c.close(); c = con.cursor()
      c.execute(query,(SNN,))
      snn_gloc = c.fetchone()
      c.close()
      # print VAR,GNN,SNN,var_gloc,snn_gloc
      if var_gloc[0] != snn_gloc[0]:
        gd2snn   = 'NA'
      else:
        gd2snn   = abs(var_gloc[1] - snn_gloc[1])
    # Calculate SD to GNN
    c = con.cursor()
    query  = "SELECT a.structid,biounit,model,x,y,z "
    query += "FROM Residue as a "
    query += "INNER JOIN GenomicIntersection as b ON a.structid=b.structid "
    query += "AND a.chain=b.chain AND a.seqid=b.seqid "
    query += "INNER JOIN GenomicConsequence as c ON b.gc_id=c.gc_id "
    query += "WHERE c.name=%s "
    c.execute(query,(VAR,))
    var_slocs = [list(i) for i in c.fetchall()]
    c.close(); c = con.cursor()
    c.execute(query,(GNN,))
    gnn_slocs = [list(i) for i in c.fetchall()]
    sd2gnn = 'NA'
    gnn_nsloc = ['NA','NA','NA','NA','NA','NA']
    for var_sloc in var_slocs:
      for gnn_sloc in gnn_slocs:
        if var_sloc[0:2] == gnn_sloc[0:2]:
          var_sloc[2:5] = [float(i) for i in var_sloc[2:5]]
          gnn_sloc[2:5] = [float(i) for i in gnn_sloc[2:5]]
          sdist = math.sqrt((var_sloc[3]-gnn_sloc[3])**2
                           +(var_sloc[4]-gnn_sloc[4])**2
                           +(var_sloc[5]-gnn_sloc[5])**2)
          if sdist < sd2gnn:
            sd2gnn   = sdist
            gnn_nsloc = gnn_sloc
    row.extend(list(gnn_nsloc))
    row.append(sd2gnn)    
    row.extend(list(snn_gloc))
    row.append(gd2snn)
    writer.writerow(row)
fout.close()
