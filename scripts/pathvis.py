# Chimera executable
if __name__ == '__main__':
  from chimera import runCommand as rc
  from chimera import replyobj
  import sys,os
  args = sys.argv
  params = {"pdbf":args[1],"neutral":args[2],"pathogenic":args[3],
            "variant":args[4],"neutcon":args[5],"pathcon":args[6],
            "pathprox":args[7],"out":args[8]}

  # Visualize with Chimera
  # Assign the annotation to relevant residues
  rc("open %(pdbf)s"%params)
  # Initialize the settings
  rc("background solid white")
  rc("ribbackbone")
  rc("~disp")
  rc("ribspline card smooth strand")
  rc("represent bs")
  rc("setattr m autochain 0")
  rc("setattr m ballScale .7")
  rc("color grey,r")
  rc("transparency 70,r")
  # Emphasize pathogenic if available, else neutral variants
  if os.path.exists(params['pathogenic']):
    rc("defattr %(pathogenic)s raiseTool false"%params)
    rc("define plane name p1 :/pathogenic")
  else:
    rc("defattr %(neutral)s raiseTool false"%params)
    rc("define plane name p1 :/neutral")
  rc("align p1")
  rc("~define")
  rc("center")
  # rc("scale 1.2")
  # Display neutral only
  if os.path.exists(params['neutral']):
    rc("defattr %(neutral)s raiseTool false"%params)
    rc("disp :/neutral & @ca")
    rc("color blue,a")
    rc("copy file %(out)s_neutral.png width 3 height 3 units inches dpi 300"%params)
    rc("~disp")
  # Display pathogenic only
  if os.path.exists(params['pathogenic']):
    rc("defattr %(pathogenic)s raiseTool false"%params)
    rc("disp :/pathogenic & @ca")
    rc("color red,a")
    rc("copy file %(out)s_pathogenic.png width 3 height 3 units inches dpi 300"%params)
    rc("~disp")
  # Display both pathogenic and neutral
  if os.path.exists(params['variant']):
    rc("defattr %(variant)s raiseTool false"%params)
    rc("disp :/pathogenicity>-9 & @ca")
    rc("rangecolor pathogenicity,a 1 red 0 blue")
    rc("copy file %(out)s_variants.png width 3 height 3 units inches dpi 300"%params)
    rc("~disp")
  rc("transparency 0,r") # confirm opacity
  # Display neutral constraint
  if os.path.exists(params['neutcon']):
    rc("defattr %(neutcon)s raiseTool false"%params)
    rc("rangecolor neutcon,r max blue min white")
    rc("copy file %(out)s_neutcon.png width 3 height 3 units inches dpi 300"%params)
  # Display pathogenic constraint
  if os.path.exists(params['pathcon']):
    rc("defattr %(pathcon)s raiseTool false"%params)
    rc("rangecolor pathcon,r max red min white")
    rc("copy file %(out)s_pathcon.png width 3 height 3 units inches dpi 300"%params)
  # Display PathProx
  if os.path.exists(params['pathprox']):
    rc("defattr %(pathprox)s raiseTool false"%params)
    rc("rangecolor pathprox,r max red min blue 0 white")
    rc("copy file %(out)s_pathprox.png width 3 height 3 units inches dpi 300"%params)
  # Export the scene and exit
  rc("disp :/pathogenicity>-9 & @ca")
  rc("rangecolor pathogenicity,a 1 red 0 blue")
  rc("save ./%(out)s.py"%params)
  rc("close all")
  rc("stop now")
