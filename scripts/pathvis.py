# Chimera executable
if __name__ == '__main__':
  from chimera import runCommand as rc
  from chimera import replyobj
  import sys,os
  args = sys.argv
  label = args[1]
  sid   = args[2]
  bio   = args[3]

  pdbf     = "%s_%s_%s.pdb"%(label,sid,bio)      # pdb file
  pdbrf    = "%s_%s_%s_renum.pdb"%(label,sid,bio)# pdb file (renumbered)
  nvattrf  = "%s_neutral.attr"%label             # neutral variants
  nvrattrf = "%s_renum_neutral.attr"%label       # neutral variants (renumbered)
  pvattrf  = "%s_pathogenic.attr"%label          # pathogenic variants
  pvrattrf = "%s_renum_pathogenic.attr"%label    # pathogenic variants (renumbered)
  avattrf  = "%s_variants.attr"%label            # all variants
  avrattrf = "%s_renum_variants.attr"%label      # all variants (renumbered)
  ncattrf  = "%s_neutcon.attr"%label             # neutral constraint
  ncrattrf = "%s_renum_neutcon.attr"%label       # neutral constraint (renumbered)
  pcattrf  = "%s_pathcon.attr"%label             # pathogenic constraint
  pcrattrf = "%s_renum_pathcon.attr"%label       # pathogenic constraint (renumbered)
  ppattrf  = "%s_pathprox.attr"%label            # path prox
  pprattrf = "%s_renum_pathprox.attr"%label      # path prox (renumbered)
  out      = "%s"%label                          # output label
  outr     = "%s_renum"%label                    # output label (renumbered)

  def visualize(pdbf,nvattrf,pvattrf,avattrf,ncattrf,pcattrf,ppattrf,out,png=True):
    # Visualize with Chimera
    # Assign the annotation to relevant residues
    rc("open %s"%pdbf)
    # Initialize the settings
    rc("background solid white")
    rc("ribbackbone")
    rc("~disp")
    rc("ribspline card smooth strand")
    rc("represent bs")
    rc("setattr m autochain 0")
    rc("setattr m ballScale .7")
    rc("color grey,r")
    # Align to neutral variants so that ClinVar and COSMIC match
    # If neutral variants unavailable, use the default orientation
    if os.path.exists(nvattrf):
      rc("defattr %s raiseTool false"%nvattrf)
      rc("define plane name p1 :/neutral")
    rc("align p1")
    rc("~define")
    rc("center")
    # rc("window")
    # Display structure only
    if png:
      rc("copy file %s_structure.png width 3 height 3 units inches dpi 300"%out)
    # Reduce ribbon transparency for variant plots
    rc("transparency 70,r")
    # Display neutral only
    if os.path.exists(nvattrf) and png:
      rc("defattr %s raiseTool false"%nvattrf)
      rc("disp :/neutral & @ca")
      rc("color blue,a")
      rc("copy file %s_neutral.png width 3 height 3 units inches dpi 300"%out)
      rc("~disp")
    # Display pathogenic only
    if os.path.exists(pvattrf) and png:
      rc("defattr %s raiseTool false"%pvattrf)
      rc("disp :/pathogenic & @ca")
      rc("color red,a")
      rc("copy file %s_pathogenic.png width 3 height 3 units inches dpi 300"%out)
      rc("~disp")
    # Display both pathogenic and neutral
    if os.path.exists(avattrf) and png:
      rc("defattr %s raiseTool false"%avattrf)
      rc("disp :/pathogenicity>-9 & @ca")
      rc("rangecolor pathogenicity,a 1 red 0 blue")
      rc("copy file %s_variants.png width 3 height 3 units inches dpi 300"%out)
      rc("~disp")
    rc("transparency 0,r") # confirm opacity
    # Display neutral constraint
    if os.path.exists(ncattrf) and png:
      rc("defattr %s raiseTool false"%ncattrf)
      rc("rangecolor neutcon,r min red max white")
      rc("copy file %s_neutcon.png width 3 height 3 units inches dpi 300"%out)
    # Display pathogenic constraint
    if os.path.exists(pcattrf) and png:
      rc("defattr %s raiseTool false"%pcattrf)
      rc("rangecolor pathcon,r max red min white")
      rc("copy file %s_pathcon.png width 3 height 3 units inches dpi 300"%out)
    # Display PathProx
    if os.path.exists(ppattrf) and png:
      rc("defattr %s raiseTool false"%ppattrf)
      rc("rangecolor pathprox,r max red min blue 0 white")
      rc("copy file %s_pathprox.png width 3 height 3 units inches dpi 300"%out)
    # Export the scene and exit
    rc("disp :/pathogenicity>-9 & @ca")
    rc("rangecolor pathogenicity,a 1 red 0 blue")
    rc("save ./%s.py"%out)
    os.chmod("./%s.py"%out,774)
    rc("close all")

  # Visualize with original PDB numbering
  print "\nChimera visualization using original PDB numbering..."
  visualize(pdbf,nvattrf,pvattrf,avattrf,ncattrf,pcattrf,ppattrf,out)

  # Visualize with renumbered PDB
  print "\nChimera visualization using renumbered PDB..."
  visualize(pdbrf,nvrattrf,pvrattrf,avrattrf,ncrattrf,pcrattrf,pprattrf,outr,png=False)

  rc("stop now")