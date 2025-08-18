test_watermark <- function(imgpath, text="百奥智汇", color="red", out_path, fig_type="png"){
  library(magick)
  
  imgpath <- normalizePath(imgpath)
  img_filename <- basename(imgpath)
  out <- file.path(out_path, img_filename)
  
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }
  
  img <- image_read(imgpath)
  img_info <- image_info(img)
  img_w <- img_info$width
  img_h <- img_info$height
  
  cat("Image:", img_w, "x", img_h, "\n")
  
  # Large, visible watermark
  size <- round(min(img_w, img_h)/8)
  cat("Size:", size, "\n")
  
  # Simple test - one watermark in center
  watermarked_image <- img %>% 
    image_annotate(
      text = text,
      color = color,
      size = size,
      font = "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
      location = paste0("+", round(img_w/2), "+", round(img_h/2)),
      degrees = 0
    )
  
  image_write(watermarked_image, out, format = fig_type)
  cat("Saved to:", out, "\n")
}