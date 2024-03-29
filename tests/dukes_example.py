# dukes_example.py

import sys
import dukes

if __name__ == '__main__':

    #Tx,mx = 5,1
    Tx = float(sys.argv[1])  # DM kinetic energy, MeV
    mx = float(sys.argv[2])  # DM mass, MeV
    vx = dukes.vBDM(Tx,mx)   # BDM velocity
    
    print(vx)                # Print the BDM velocity