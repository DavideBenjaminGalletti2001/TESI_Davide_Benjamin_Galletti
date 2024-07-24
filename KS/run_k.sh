
j="0.0"

while [[ "$(bc <<< "$j < 0.049")" == "1" ]]; do
    # do something with i
    j="$(bc <<< "$j + 0.001")"
    sed -i "5c ${j}" input.dat
    ./main.exe 
    python3 main.py
done