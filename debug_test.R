debug_test <- function(imgpath, out_path){
  library(magick)
  
  img <- image_read(imgpath)
  img_info <- image_info(img)
  img_w <- img_info$width
  img_h <- img_info$height
  
  cat("Original image size:", img_w, "x", img_h, "\n")
  
  # Test 1: English text without font
  test1 <- img %>% 
    image_annotate(
      text = "TEST1",
      color = "red",
      size = 100,
      location = "+100+100"
    )
  image_write(test1, file.path(out_path, "test1_english.png"))
  
  # Test 2: Chinese text without font
  test2 <- img %>% 
    image_annotate(
      text = "百奥智汇",
      color = "red", 
      size = 100,
      location = "+100+200"
    )
  image_write(test2, file.path(out_path, "test2_chinese_no_font.png"))
  
  # Test 3: Chinese with font path
  test3 <- img %>% 
    image_annotate(
      text = "百奥智汇",
      color = "red",
      size = 100,
      font = "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
      location = "+100+300"
    )
  image_write(test3, file.path(out_path, "test3_chinese_with_font.png"))
  
  cat("Created 3 test images\n")
}