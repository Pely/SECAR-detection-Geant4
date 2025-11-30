#!/bin/bash

flag=0
startRun=0
stopRun=27

#Gia na tre3ei grafeis ./run.sh #(apo pou 3ekinaei) #(pou teleiwnei)

i=${startRun}
j=100

while [[ ${flag} == 0 ]]
do


    ################################################################################
    while [[ ${i} -le ${stopRun} ]]
    do
        
        #--------------------------------------------------------------------------#
        #loop   
        for ((i = 0; i <= stopRun; i++)) 
        do  
          
          if [[ $i==0 ]];
          then
            make
          fi
          
          ./SECAR ../neutronInput/neutrons/n_${j}.mac

          hadd -f neutrons/0deg/n_${j}keV.root output_t*

          # If this is the last run then break
          if [[ ${i} == ${stopRun} ]]; 
          then
            flag=1
            break
          fi
        
          # Increment run number
          if [[ j -ge 100 && j -lt 1500 ]]; then
            let j=j+100
          elif [[ j -ge 1500 && j -lt 6000 ]]; then
            let j=j+500
          elif [[ j -ge 6000 ]]; then
            let j=j+1000
          fi

        done
      
        # This is somewhat redundant, but it works.
        if [[ ${flag} == 1 ]]; 
        then
          break
        fi

        #--------------------------------------------------------------------------#
        
        wait
    done
    ##################################################################################

    wait
done

