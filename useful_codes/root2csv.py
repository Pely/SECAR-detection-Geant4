import ROOT
import csv
import math

# --- SETTINGS ---
root_file = "bld_an/ne-20_pg.root"
tree_name = "Target"         # change to your TTree name
branches_to_save = ["posX_str", "posY_str", "theta_str", "phi_str", "ekin_str"]   # choose the branches you want
csv_out = "bld_an/stripper.csv"

transforms = {
    "posX": lambda v: v * 1e-3,     # mm → m
    "posY": lambda v: v * 1e-3,     # mm → m
    # add more if needed, leave the rest untouched
}

derived_vars = {
    "ax": lambda ev: math.tan(ev.theta) * math.cos(ev.phi),
    "ay": lambda ev: math.tan(ev.theta) * math.sin(ev.phi),
}

# --- OPEN FILE & TREE ---
f = ROOT.TFile.Open(root_file)
t = f.Get(tree_name)

# --- WRITE CSV ---
with open(csv_out, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    
    # Write header row
    #writer.writerow(branches_to_save)
    writer.writerow(["x_m", "y_m", "theta_rad", "phi_rad", "e_MeV", "ax_rad","ay_rad"])
    
    # Loop over events
    for event in t:
        row = []
        for br in branches_to_save:
            val = getattr(event, br)

            # apply transform if needed
            if br in transforms:
                val = transforms[br](val)
            row.append(val)
            
        for name, func in derived_vars.items():
            row.append(func(event))
            
        writer.writerow(row)


print("Saved CSV:", csv_out)
