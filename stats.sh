#!/bin/bash

# Copyright 2023 Daniel Grier and Luke Schaeffer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

if [[ -n "$1" ]]; then
    let step=$1
else
    let step=2
fi

printf '' > stats.csv

for ((j=$step; ; j+=$step))
do
    if [[ -n $2 ]] && [[ $j -gt $2 ]]; then
        break
    fi
    printf "Timing %dx%d grid... " $j $j

    printf "%d," $j >> stats.csv
    ./gridCHP++ -l $j >> stats.csv
    printf "\n" >> stats.csv

    printf "done!\n"
done
