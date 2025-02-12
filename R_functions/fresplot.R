fresplot <- function(fresiduals, 
                     x,
                     title=" ",
                     scale=c("uniform","normal"),
                     xl = -3,
                     xp = 3,
                     yl = NULL,
                     yp = NULL,
                     xlabs="X",
                     heatmapcut=101,
                     is.binary = FALSE) {
  scale <- match.arg(scale)
  length <- length(x) 
  # Initialize the matrix to store numbers
  numbers <- matrix(NA, nrow = length, ncol = heatmapcut)
  
  # Fill the matrix with the computed values based on fresiduals and x
  for (h in 1:length) {
    for (a in 1:ncol(numbers)) {
      numbers[h, a] <- fresiduals[h, 1] + (fresiduals[h, 2] - fresiduals[h, 1]) / (heatmapcut-1) * (a - 1)
    }
  }
  
  # Replicate x to match the number of rows in numbers
  x_101 <- rep(x, heatmapcut)
  if (is.binary==TRUE){
    x_101<-x_101+runif(length(x_101),min=0,max=0.01)
  }
  # Convert the matrix to a vector and then apply the normal quantile transformation
  numbers_v <- as.vector(numbers)
  qnumbers_v <- qnorm(numbers_v)
  
  # Combine the x11 and transformed values into a data frame
  qnumbers <- cbind.data.frame(x_101, qnumbers_v)
  numbers <- cbind.data.frame(x_101, numbers_v)
  
  if (scale=="uniform"){
    # Create the plot for the uniform scale (functional residuals)
    yl <- ifelse(is.null(yl), 0, yl)
    yp <- ifelse(is.null(yp), 1, yp)
    p_unif <- ggplot(numbers, aes(x_101, numbers_v)) + 
      stat_density_2d(aes(fill = stat(level)), geom = 'polygon') + 
      scale_fill_viridis_c(name = "density") + 
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
      geom_smooth(method = "loess", se = FALSE) + 
      xlim(xl,xp) +  # Set x-axis limits
      ylim(yl,yp) +
      xlab("X") + 
      ylab("") + 
      labs(title = title) + 
      theme(plot.title = element_text(size = 12), 
            axis.title = element_text(size = 12))
    # Return the plot
    return(p_unif)
  }else if (scale == "normal") {
    # Normal scale plot
    if(is.null(yl) & is.null(yp)){
    p_norm<-ggplot(qnumbers, aes(x_101,qnumbers_v)) +
      stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
      scale_fill_viridis_c(name = "density")+
      geom_hline(yintercept=0,linetype="dashed", color = "red")+
      labs(x = xlabs, y=" ")+
      xlim(xl,xp) + # Set x-axis limits 
      geom_smooth(method = "loess",se=FALSE)+
      labs(title = title)+
      theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
    }else {
      p_norm<-ggplot(qnumbers, aes(x_101,qnumbers_v)) +
        stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
        scale_fill_viridis_c(name = "density")+
        geom_hline(yintercept=0,linetype="dashed", color = "red")+
        labs(x = xlabs, y=" ")+
        xlim(xl,xp) + ylim(yl, yp) + # Set x-axis limits 
        geom_smooth(method = "loess",se=FALSE)+
        labs(title = title)+
        theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
    }
    # Return the plot
    return(p_norm)
  }}
