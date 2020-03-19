
volcano_plot <- function (dat, out_dir, out_name, title, fc_line, alpha){
  res_df <- as.data.frame(dat)
  res_df$sig <- c()
  res_df$cols <- c()
  for(i in 1:nrow(res_df)){
    if(is.na(res_df$padj[i])){
      res_df$sig[i] <- NA
      res_df$cols[i] <- NA
    }
    else if(res_df$padj[i]<=alpha & res_df$log2FoldChange[i] > fc_line){
      res_df$sig[i] <- paste0("Q-value <", alpha, ", Log2FC >", fc_line)
      res_df$cols[i] <- "indianred"
    } 
    else if(res_df$padj[i]<=alpha & res_df$log2FoldChange[i] < -fc_line){
      res_df$sig[i] <- paste0("Q-value <", alpha, ", Log2FC < -", fc_line)
      res_df$cols[i] <- "indianred"
    } 
    else if(res_df$padj[i]<=alpha & res_df$log2FoldChange[i]>-fc_line & res_df$log2FoldChange[i]<fc_line){
      res_df$sig[i] <- paste0("Q-value <", alpha, ", Log2FC > -", fc_line, "and ", ", Log2FC <", fc_line)
      res_df$cols[i] <- "cornflowerblue"
    } 
    else if(res$padj[i]>alpha & res_df$log2FoldChange[i] > fc_line){
      res_df$sig[i] <- paste0("Q-value >", alpha, ", Log2FC >", fc_line)
      res_df$cols[i] <- "gray47" 
    }
    else if(res$padj[i]>alpha & res_df$log2FoldChange[i] < -fc_line){
      res_df$sig[i] <- paste0("Q-value >", alpha, ", Log2FC < -", fc_line)
      res_df$cols[i] <- "gray47" 
    }
    else if(res_df$padj[i]>alpha & res_df$log2FoldChange[i] < fc_line){
      res_df$sig[i] <- paste0("Q-value >", alpha, ", Log2FC < -", fc_line)
      res_df$cols[i] <- "gray10" 
    }
  }
  
  ppi=300
  png(paste0(out_dir, out_name), width=6*ppi, height=7*ppi, res=ppi)
  p = ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(col=sig),alpha = 0.5, shape = 21, size =2.5, colour = res_df$cols, fill = res_df$cols)  + 
    #xlim(-2.5, 2.5) +
    xlab("Log2 fold change estimate") + ylab("-log10(Q-value)") +
    geom_hline(yintercept = -log10(alpha), color = "black", linetype = "dashed", size = 0.4) + 
    labs(color="") + 
    theme(legend.key = element_blank()) + 
    ggtitle(title) + 
    geom_vline(xintercept = fc_line, colour = "black", linetype="dotted") + 
    geom_vline(xintercept = -fc_line, colour = "black", linetype="dotted") + 
    theme(legend.key.size = unit(1, "cm"), 
          #panel.grid.major = element_line(colour = "#d3d3d3"), 
          panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
          panel.background = element_blank(),
          axis.text.x=element_text(colour="black", size = 20, hjust = 1), 
          axis.text.y=element_text(colour="black", size = 20), 
          axis.title.x=element_text(colour="black", size = 20), 
          axis.title.y=element_text(colour="black", size = 20), 
          legend.key = element_blank(), 
          legend.text = element_text(colour="black", size = 20), 
          legend.title = element_text(colour="black", size = 20))
  print(p)
  dev.off()
  
}
