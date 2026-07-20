# campaign-finance-networks

In this project I model FEC campaign finance data as an absorbing Markov chain in which candidates are the absorbing states. As such, to "follow the money" is to follow the random walk of money. But most committees' political affiliation are unknown. Using this model, I (1) empirically infer entity ideologies and subsequently (2) compare Democratic and Republican random walk dynamics.  

My takeaways: Individual ideology is an emergent feature of the whole network. While Democratic and Republican dollars take paths of similar lengths, Republican paths end in more in-party candidates (despite Democrats having more candidates overall) and more easily cross the political boundary than Democratic ones. 

To run my data pipeline, run `data/download_data.sh` (this downloads and cleans the 2000-2022 FEC data) then `data/process_data.sh` (this prepares the data for use by `igraph`). Run the analyses in `analysis/analysis.ipynb`. Note that the null model simulations use `run_null_simulations.py` which takes around 45 minutes to run. 

<img width="1019" height="972" alt="absorbing" src="https://github.com/user-attachments/assets/6fd0222b-d74e-49e3-9218-a884ab95c5fd" />

This project was done at the Complexity Science Hub in Vienna as part of their summer internship program, which I found through Princeton. 
