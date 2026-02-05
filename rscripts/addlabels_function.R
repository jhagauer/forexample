# Function to add ecotype, treatment and replicate information
# Might rewrite to make use of sample_info.tsv (also better if information had to be added or changed)
add_labels <- function(x){

        HSdown <-c("sample01", "sample06","sample13", "sample16", "sample18", "sample20")
        HSup <-c("sample03", "sample05", "sample09", "sample11",  "sample21", "sample23")
        LSdown <-c("sample04", "sample07", "sample10", "sample14",  "sample17", "sample24")
        LSup <-c("sample02", "sample08", "sample12", "sample15","sample19", "sample22")


        repB <- c("sample09","sample10","sample11","sample12","sample13","sample14","sample15","sample16")
        repA <- c("sample01","sample02","sample03","sample04","sample05","sample06","sample07","sample08")
        repC <- c("sample21","sample22","sample23","sample24","sample20","sample19","sample18","sample17")

        down <- c(HSdown, LSdown)
        up <- c(HSup, LSup)

        LS <- c(LSup, LSdown)
        HS <- c(HSup, HSdown)

        x[["bio_rep"]] <- x[["orig.ident"]]

        for(sample in repA){
            x[["bio_rep"]][which(x[["bio_rep"]] == sample),] <- "A"
        }
        for(sample in repB){
            x[["bio_rep"]][which(x[["bio_rep"]] == sample),] <- "B"
        }
        for(sample in repC){
            x[["bio_rep"]][which(x[["bio_rep"]] == sample),] <- "C"
        }

        x[["group"]] <- x[["orig.ident"]]

        for(sample in HSdown){
            x[["group"]][which(x[["group"]] == sample),] <- "HSdown"
        }
        for(sample in HSup){
            x[["group"]][which(x[["group"]] == sample),] <- "HSup"
        }
        for(sample in LSdown){
            x[["group"]][which(x[["group"]] == sample),] <- "LSdown"
        }
        for(sample in LSup){
            x[["group"]][which(x[["group"]] == sample),] <- "LSup"
        }

        x[["treat"]] <- x[["orig.ident"]]

        for(sample in LS){
            x[["treat"]][which(x[["treat"]] == sample),] <- "LS"
        }
        for(sample in HS){
            x[["treat"]][which(x[["treat"]] == sample),] <- "HS"
        }

        x[["ecotype"]] <- x[["orig.ident"]]

        for(sample in up){
            x[["ecotype"]][which(x[["ecotype"]] == sample),] <- "up"
        }
        for(sample in down){
            x[["ecotype"]][which(x[["ecotype"]] == sample),] <- "down"
        }

        return(x)
}