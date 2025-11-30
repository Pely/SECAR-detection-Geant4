import ROOT
import pickle

# === USER SETTINGS ===
input_file  = "bld_an/ne-20_5107_bigfoil.root"
tree_name   = "Target"
output_file = "na-21_rays_afterCfoil.pkl"

# Branch names in your ROOT tree:
branch_x = "posX"
branch_y = "posY"
branch_theta = "theta"   # polar angle
branch_phi   = "phi"     # azimuthal angle
branch_E     = "ekin"       # energy (not used now)
branch_Z = "Z"
branch_A = "A"

# ======================

# Open file and tree
f = ROOT.TFile.Open(input_file)
t = f.Get(tree_name)

rays = []

for i in range(t.GetEntries()):
    t.GetEntry(i)
    a = int(getattr(t, branch_A))
    z = int(getattr(t, branch_Z))
    if(z==11 and a==21):
        x = float(getattr(t, branch_x))
        y = float(getattr(t, branch_y))
        theta = float(getattr(t, branch_theta))
        phi = float(getattr(t, branch_phi))
        e = float(getattr(t, branch_E))

        # === Convert theta, phi → ax, by ===
        # ax = tan(theta) * cos(phi)
        # by = tan(theta) * sin(phi)
        # (Assuming small angles, ax ≈ theta*cos(phi), by ≈ theta*sin(phi))
        ax = theta * ROOT.TMath.Cos(phi)
        by = theta * ROOT.TMath.Sin(phi)
        dE = (e - 21.7)/21.7
        l = dm = dz = 0.0

        rays.append(([x], [ax], [y], [by], [l], [dE], [dm], [dz]))
        print(rays)
# Save to pickle
with open(output_file, "wb") as f_out:
    pickle.dump(rays, f_out, protocol=pickle.HIGHEST_PROTOCOL)

print(f"Saved {len(rays)} rays to {output_file}")
