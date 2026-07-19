# campaign-finance-networks

In this project I model campaign finance as an absorbing Markov chain. As such, "follow the money" means "follow the random walk of money." Using this model, I empirically infer entity ideologies and subsequently compare Democratic and Republican random walk dynamics.  

To run my data pipeline, run `data/download_data.sh` (this downloads and cleans the 2000-2022 FEC data) then `data/process_data.sh` (this prepares the data for use by `igraph`). Run the analyses in `analysis/analysis.ipynb`. 

<img width="1019" height="972" alt="absorbing" src="https://github.com/user-attachments/assets/6fd0222b-d74e-49e3-9218-a884ab95c5fd" />
