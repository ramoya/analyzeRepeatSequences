# analyzeRepeatSequences


Input file formats:

        unitorders_df
        
        A parsed multiple sequence alignment. Each row is one repeat unit. Rows are ordered by individual.
        Within an individual, rows are ordered from bottom to top in the 5' to 3' direction.
        
        Columns:
              character: repeat unit encoded as a single ASCII character (character)
              sample: individual sequence ID (character)
              count: width of each repeat unit on the plot (numeric)
              seq: sequence of repeat unit
              gp: integer that increments one every row over the ungrouped dataframe. 
                  This gives every repeat unit in a sequence its own position, 
                  rather than clustering identical units together (integer)
              group: used with plotSubgroups = T to cluster sequences in the plot
              
To install:

        devtools::install_github("ramoya/analyzeRepeatSequences")
              
              
Preprint: https://www.medrxiv.org/content/10.1101/2024.03.05.24303780v2

              

