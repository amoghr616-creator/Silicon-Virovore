load "/Users/centurion616/Desktop/Silicon-Virovore/prototype4/data/haddock3/receptors/9sya_source_chain_A.pdb", receptor
wizard mutagenesis
get_wizard().set_mode("CYS")
get_wizard().do_select("chain A and resi 141")
get_wizard().apply()
get_wizard().set_mode("THR")
get_wizard().do_select("chain A and resi 167")
get_wizard().apply()
get_wizard().set_mode("ILE")
get_wizard().do_select("chain A and resi 185")
get_wizard().apply()
python
from pymol import cmd
cmd.set_wizard()
cmd.save("/Users/centurion616/Desktop/Silicon-Virovore/prototype4/data/haddock3/receptors/9sya_reference_corrected_chain_A.pdb", "receptor")
python end
quit
