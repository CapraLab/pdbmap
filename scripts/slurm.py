#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : slurm.py
# Authors        : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-01-22
# Description    : Launches and tracks SLURM submissions.
#=============================================================================#

import subprocess as sp
import numpy as np
import sys,os
from time import sleep
from copy import deepcopy

class SlurmJob():

  def __init__(self,script,args=[],array=None,account=None,debug=False):
    self.script   = script  # File name or content of slurm script
    if not os.path.exists(script):
      self.script = self._create_temp(self.script)
    self.args     = args    # Arguments to the slurm script
    self.array    = array   # Array string (e.g. "0-99")
    self.account  = account # SLURM account name (if not user default)
    self.debug    = debug   # If true, submit job to debug partition
    # Submit the job
    self.jobid    = self._submit()
    self.failed   = False # Fail flag for job submission failures
    self.finished = False
    self.info     = None

  def _create_temp(self,script):
    """ Writes the SLURM contents to a temporary file for submission """
    fname = "../temp/psb_pathprox_%d.slurm"%np.random.uniform(100000,999999)
    with open(fname,'wb') as fout:
      fout.write(script)
    return fname

  def _submit(self):
    """ Submits the  job to SLURM"""
    job_submit = ["sbatch"]
    if self.account:
      job_submit.extend(["--account",self.account])
    if self.array:
      if self.debug:
        job_submit.extend(["--array","0-1"])
      else:
        job_submit.extend(["--array",self.array])
    if self.debug:
      job_submit.extend(["--time","30","--partition","debug"])
    job_submit.append(self.script)
    job_submit.extend(self.args)
    p = sp.Popen(job_submit,stdout=sp.PIPE,stderr=sp.PIPE)
    # Extract the job ID
    stdout,stderr = p.communicate()
    # Extract and return the jobid (e.g. "Submitted batch job 12345678")
    if stderr:
      raise Exception("%s\n"%stderr)
    return stdout.strip().split()[-1]

  def get_info(self):
    if self.info:
      return self.info
    else:
      show_job = ["scontrol","show","job",self.jobid]
      p = sp.Popen(show_job,stdout=sp.PIPE,stderr=sp.PIPE)
      stdout,stderr = p.communicate()
      return dict(tuple(info.split('=')) for info in stdout.split() if len(info.split('='))==2)

  def get_state(self):
    if self.failed:
      return "FAILED"
    return self.get_info()["JobState"]

  def get_exit_status(self):
    return self.get_info()["ExitCode"]

  def get_stdout(self):
    with open(self.get_info()["StdOut"],'rb') as fin:
      return fin.readlines()

  def get_stderr(self):
    with open(self.get_info()["StdErr"],'rb') as fin:
      return fin.readlines()

  def get_submit_time(self):
    return self.get_info()["SubmitTime"]

  def get_start_time(self):
    return self.get_info()["StartTime"]

  def get_end_time(self):
    return self.get_info()["EndTime"]

  def is_pending(self):
    return self.get_state() == "PENDING"

  def is_running(self):
    return self.get_state() == "RUNNING"

  def is_complete(self):
    return self.get_state() == "COMPLETED"

  def is_failed(self):
    return self.get_state() == "FAILED"

  def is_cancelled(self):
    return self.get_state() == "CANCELLED"

  def is_timeout(self):
    return self.get_state() == "TIMEOUT"

  def is_finished(self):
    if not self.finished:
      self.finished = self.get_state() in \
        ["TIMEOUT","SPECIAL_EXIT","PREEMPTED","NODE_FAIL",
         "FAILED","COMPLETED","CANCELLED","BOOT_FAIL"]
      if self.finished:
        # Record job exit info and stop polling SLURM
        self.info = self.get_info()
    return self.finished