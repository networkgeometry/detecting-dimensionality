# Script to generate SD surrogates of a given network. 
# 
# The script requires an edgelist of the given network and a parameter setting.
# 
# Parameters:
# 
# - network: name of the network (edgelist file must be located in RealNets folder with edge extension and features file in RealFeats with csv extension)
# - resolution: number of surrogates per dimension
# - n_poll: number of points used to infer the relation T vs. Beta
# - wsize: size of the clustering interval in which to create the surrogates
# - nrealizations: number of realizations per random Beta value (default=1)  
# - maxD: the script will generate surrogates from D=1 to D=maxD
# 
# The resulting surrogates will be placed in SDnets folder. Some folders will be created during the process for calculation purposes.
# 
########### PART 1: Infering relation T (proportion of triangles) vs. Beta ###########

# Creating necessary folders

mkdir RealFeats
mkdir SDnets
mkdir SDnets/$1
mkdir temp
mkdir errors
mkdir outputs
# Recieving parameters

network=$1
RES=$2
step_sondeo=$3
wsize=$4
nrealizations=$5
maxD=$6

# Calculating T, maxT (T + wsize/2) and minT (T + wsize/2)

d=1
beta=$(echo " $d + 0.5 " |bc -l)
min=1

Bstep=$(echo " 13.5 / $step_sondeo " |bc -l)

iter=1

FILE="RealFeats/"$network".csv" 
if ! test -f "$FILE"; then
    ./task_cyclesmap.sh "RealNets/"$1".edge" $FILE
fi
REALT=$( awk -v line=1 'FNR==line{print $1}'  $FILE )
minT=$(echo "$REALT - ($wsize / 2)" |bc -l) 
maxT=$(echo "$REALT + ($wsize / 2)" |bc -l) 
echo "Proportion of triangles in netowrk: "$REALT
echo "Minimum value of T for surrogates: "$minT
echo "Maximum value of T for surrogates: "$maxT

# Declaring a dictionary betaT to store relaions T vs. Beta Beta-Clustering Dictionary for D=1 

declare -A betaT

# Generating n_poll equidistant surrogates in dim=1 to infer the relation T vs. B in dim=1 

./task_SD.sh "RealNets/"$1 "temp/"$1"_"$d"_"$beta"_1" $d $beta 1

while (( $(echo "$beta < 15" |bc -l) )); #Loop over beta values (using a x100 factor) generating instances
do 
   beta=$(echo " $beta + $Bstep " |bc -l) 
   ./task_SD.sh "RealNets/"$1 "temp/"$1"_"$d"_"$beta"_1" $d $beta 1
done

# Calculating proportion of triangles (T) of those n_poll surrogates

beta=$(echo " $d + 0.5 " |bc -l)
FILE="temp/"$1"_"$d"_"$beta"_1.csv" 
./task_triangles.sh "temp/"$1"_"$d"_"$beta"_1.edge" $FILE
while (( $(echo "$beta < 15" |bc -l) )); 
do 
   beta=$(echo " $beta + $Bstep " |bc -l) 
   FILE="temp/"$1"_"$d"_"$beta"_1.csv" 
   ./task_triangles.sh "temp/"$1"_"$d"_"$beta"_1.edge" $FILE
done

# Storing results and printing

echo "Beta-Clustering Dictionary for Dim 1"

beta=$(echo " $d + 0.5 " |bc -l)
while (( $(echo "$beta < 15" |bc -l) )); 
do 
   FILE="temp/"$1"_"$d"_"$beta"_1.csv" 
   T=$( awk -v line=1 'FNR==line{print $1}'  $FILE )
   betaT[$beta]=$T
   echo "Beta: "$beta
   echo "T: "${betaT[$beta]}
   beta=$(echo " $beta + $Bstep " |bc -l)
done

# Formatting T's and Beta's in order to fit the curve with python script

betas=""
clusterings=""
flag=0
for i in "${!betaT[@]}"
do
    zero=0;
    if [ $flag -eq $zero ]; then
        betas=$i
        clusterings=${betaT[$i]}
        flag=1
    else
        betas=$betas","$i
        clusterings=$clusterings","${betaT[$i]}
    fi
done

echo $betas
echo $clusterings

# Fitting curve T vs. Beta (for dim=1) with python script

outputString=($(python curve_fitting.py $betas $clusterings | tr -d '[],'))

# Printing obtained parameters of the polinomial curve

echo "c_inf:"${outputString[0]}
echo "a:"${outputString[1]}
echo "beta_min:"${outputString[2]}
echo "Params:"${outputString[@]}

########### PART 2: GENERATING SURROGATES ###########

# Flag to indicate when we reach maximum D (T exceeds maximum value)

flag=1

# Iterating over D

for d in $(seq 1 $maxD); 
do 
    # Calculating maximum T
    maxB=$(echo " 15 * $d " |bc -l) 
    if ! [ "$flag" -eq "0" ]; then
        ./task_SD.sh "RealNets/"$1 "temp/"$1"_"$d"_"$maxB"_1" $d $maxB 1
        FILE="temp/"$1"_"$d"_"$maxB"_1.csv"
        ./task_triangles.sh "temp/"$1"_"$d"_"$maxB"_1.edge" $FILE
        maxctemp=$( awk -v line=1 'FNR==line{print $1}'  $FILE )
        # Saving Maximum T for D=1
        if [[ "$d" == 1 ]]; then
            maxtemp=$( awk -v line=1 'FNR==line{print $1}'  $FILE )
        fi
        echo "maximum clustering for D="$d": "$maxctemp
        if (( $(echo "$maxctemp > $REALT" |bc -l) )); then
            # Calculating max and min values of Beta
            ratio=$(echo "($REALT * $maxtemp)/ $maxctemp" |bc -l) 
            minratio=$(echo "($minT * $maxtemp)/ $maxctemp" |bc -l) 
            maxratio=$(echo "($maxT * $maxtemp)/ $maxctemp" |bc -l) 
            echo "ratio (REALT * maxcd=1) / maxtemp = "$ratio
            echo "ratio (minT * maxcd=1) / maxtemp = "$minratio
            echo "ratio (maxT * maxcd=1) / maxtemp = "$maxratio

            mkdir SDnets/"$1"/D"$d"

            p=0
            max=$maxratio

            echo "RealT "$REALT
            echo "minratio "$minratio
            echo "maxratio "$maxratio
            if (( $(echo "$maxratio > ${outputString[0]}" |bc -l) )); then
               maxratio=${outputString[0]}
            fi            
            #Loop over resolution
            while (( $(echo "$p < $RES" |bc -l) )); 
            do
               n=1
               beta=0
               #Generating surrogates wuth random T (Beta < d + 0.05)
               while (( $(echo "$beta < ($d + 0.05)" |bc -l) )); do
                   c=$(python -c "import random;print(random.uniform("$minratio","$maxratio"))")
                   beta=$(python curve.py $c ${outputString[0]} ${outputString[1]} ${outputString[2]} )
                   echo "Beta original="$beta
                   beta=$(python -c "print("$beta"*"$d")")
               done
               echo "D="$d
               echo "Surrogate number "$p
               echo "T="$c
               echo "Beta="$beta
    
               while [ $n -le $nrealizations ] #Loop over number of realizations
               do    
                   randomico=$RANDOM
                   ./task_SD.sh "RealNets/"$1 "SDnets/"$1"/D"$d"/"$1"_"$d"_"$beta"_"$n $d $beta $randomico
                   n=$(( n + 1 ))
               done
               p=$(( p + 1 ))
            done

        else
            flag=0
        fi    
    fi
done
