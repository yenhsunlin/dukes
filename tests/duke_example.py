# duke_example.py

import sys
import duke

if __name__ == '__main__':

    #Tx,mx = 5,1
    Tx = int(sys.argv[1])   # DM kinetic energy, MeV
    mx = int(sys.argv[2])   # DM mass, MeV
    vx = duke.vBDM(Tx,mx)   # BDM velocity
    
    print(vx)               # Print the BDM velocity