ogel$set('public','log',function(...){
    cat(format(as.POSIXct(Sys.time()), format='%Y-%m-%d %H:%M:%OS '), ...)
    invisible(self)
})

ogel$set('public','sep_line',function(){
  cat('----------------------------------------\n')
})

ogel$set('public','add_watermarker',function(imgpath,text="百奥智汇",color="#696969",out_path=NULL,fig_type="png"){
  if(is.null(out_path)) {
    out_path <- self$watermark_dir
  }
  library(magick)
  library(stringr)
  library(parallel)
  
  self$say(paste("Starting watermark process for:", basename(imgpath)))
  
  # Auto-add .png extension if no extension provided
  if(!grepl("\\.", basename(imgpath))) {
    imgpath <- paste0(imgpath, ".png")
    self$say("No extension detected, added .png extension")
  }
  
  # Check if image exists, if not try fig_dir
  if(!file.exists(imgpath)) {
    self$say(paste("Image not found at:", imgpath))
    alt_path <- file.path(self$fig_dir, basename(imgpath))
    self$say(paste("Searching in fig_dir:", alt_path))
    if(file.exists(alt_path)) {
      imgpath <- alt_path
      self$say("Image found in fig_dir")
    } else {
      stop(paste("Image not found in either location:", imgpath, "or", alt_path))
    }
  } else {
    self$say("Image found at specified path")
  }
  
  imgpath <- normalizePath(imgpath)
  img_filename <- basename(imgpath)
  out <- file.path(out_path, img_filename)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }
  
  img <- image_read(imgpath)
  img_info <- image_info(img)
  img_w <- img_info$width
  img_h <- img_info$height
  self$say(paste("Image dimensions:", img_w, "x", img_h))
  
  # size大小为宽度和长度最小值的1/20 px四舍五入
  size <- round(min(img_w,img_h)/20)
  self$say(paste("Watermark font size:", size))
  # if(img_w<=1500){
  #   size <- 45
  # }else if(img_w<=3000){
  #   size <- 90
  # }else{
  #   size <- 120
  # }
  
  sudoku1_w <- img_w/3/2
  sudoku1_h <- img_h/3/2
  sudoku2_w <- img_w/2
  sudoku4_h <- img_h/2
  sudoku3_w <- img_w-img_w/3/2
  sudoku7_h <- img_h-img_h/3/2
  
  # 九宫格1的中间
  watermarked_image <- img %>% 
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku1_w}+{sudoku1_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    ) %>% 
    # 九宫格2的中间
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku2_w}+{sudoku1_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    ) %>% 
    # 九宫格3的中间
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku3_w}+{sudoku1_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    ) %>% 
    # 九宫格4的中间
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku1_w}+{sudoku4_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    ) %>% 
    # 九宫格5的中间
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku2_w}+{sudoku4_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    ) %>% 
    # 九宫格6的中间
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku3_w}+{sudoku4_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    ) %>% 
    # 九宫格7的中间
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku1_w}+{sudoku7_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    ) %>% 
    # 九宫格8的中间
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku2_w}+{sudoku7_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    ) %>% 
    # 九宫格9的中间
    image_annotate(
      text = text,       # 替换为你的水印内容
      color = color,            # 水印文字颜色
      size = size,                  # 字体大小
      font = "WenQuanYi Zen Hei",
      location = str_glue("+{sudoku3_w}+{sudoku7_h}"),
      degrees = -45,
      weight = 700                # 字体粗细（可选）
    )
  self$say("Applying 9-grid watermarks...")
  image_write(watermarked_image, out, format = fig_type)
  self$say(paste("Watermarked image saved to:", out))
  invisible(self)
})

ogel$set('public','save_fig',function(fig, file_path, width = 5, height = 5, dpi = 300, fontfamily = NULL, out_dir = NULL, auto_watermark = TRUE){
  if(is.null(out_dir)) {
    out_dir <- self$fig_dir
  }
  # Use out_dir for the file path
  full_path <- file.path(out_dir, basename(file_path))
  dev_off_safe <- function() {
    devs <- dev.list()
    if (length(devs) > 0) dev.off(devs[length(devs)])
  }
  
  print_with_font <- function(fig, fontfamily) {
    if (!is.null(fontfamily)) {
      library(grid)
      pushViewport(viewport(gp = gpar(fontfamily = fontfamily)))
      if (inherits(fig, "Heatmap")) {
        draw(fig, newpage = FALSE)
      } else {
        print(fig)
      }
      popViewport()
    } else {
      if (inherits(fig, "Heatmap")) {
        draw(fig)
      } else {
        print(fig)
      }
    }
  }
  
  # Save as PDF
  self$say("Generating PDF...")
  pdf(paste0(full_path,'.pdf'), width = width, height = height)
  print_with_font(fig, fontfamily)
  dev_off_safe()
  
  # Save as PNG
  self$say("Generating PNG...")
  png(paste0(full_path,'.png'), width = width, height = height, units = 'in', res = dpi)
  print_with_font(fig, fontfamily)
  dev_off_safe()
  
  self$say(paste("Figure saved as:", paste0(full_path, ".pdf"), "and", paste0(full_path, ".png")))
  
  # Auto watermark if enabled
  if(auto_watermark) {
    self$say("Auto watermark enabled - applying watermarks...")
    
    # Watermark both PNG and PDF if they exist
    png_path <- paste0(full_path, ".png")
    pdf_path <- paste0(full_path, ".pdf")
    
    if(file.exists(png_path)) {
      self$add_watermarker(png_path)
    }
    if(file.exists(pdf_path)) {
      self$add_watermarker(pdf_path)
    }
    
    self$say("Auto watermarking completed")
  }
  
  invisible(self)
})
