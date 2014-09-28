import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;

public class TestColorReplace {

    public static void main(String[] args) {
        String input_file = args[0];
        String output_file = args[1];
        String conversion = args[2];

        try {
            BufferedImage img = colorImage(ImageIO.read(new File(input_file)), conversion);

            ImageIO.write(img, "png", new File(output_file + "_" + conversion + ".png"));
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    private static BufferedImage colorImage(BufferedImage image, String conversion) {
        int width = image.getWidth();
        int height = image.getHeight();
        WritableRaster raster = image.getRaster();

        for (int xx = 0; xx < width; xx++) {
            for (int yy = 0; yy < height; yy++) {
                int[] pixels = raster.getPixel(xx, yy, (int[]) null);

                System.out.println("original pixel: " + pixels[0] + " " + pixels[1] + " " + pixels[2]);

                //convert to sapphire
                if (conversion.equals("sapphire")) {
                    if (pixels[2] != 64 && pixels[2] != 0) {
                        pixels[2] += Math.min(pixels[2] + 50, 255);
                    }
                }

                //convert to silver
                if (conversion.equals("sapphire")) {
                    pixels[0] *= 2;
                    pixels[1] *= 2;
                    pixels[2] *= 2;
                }

                //convert to green
                if (conversion.equals("emerald")) {
                    if (pixels[1] != 68 && pixels[1] != 0) {
                        pixels[1] += Math.min(pixels[1] + 50, 255);
                    }
                }

                //convert to ruby
                if (conversion.equals("ruby")) {
                    if (pixels[0] != 68 && pixels[0] != 0) {
                        pixels[0] += Math.min(pixels[0] + 50, 255);
                    }
                }

                //convert to gold
                if (conversion.equals("gold")) {
                    if (pixels[0] != 68 && pixels[0] != 0) {
                        pixels[0] += Math.min(pixels[0] + 75, 255);
                        pixels[1] += Math.min(pixels[1] + 50, 255);
                    }
                }

                //convert to gold
                if (conversion.equals("bronze")) {
                    if (pixels[0] != 68) {
                        pixels[0] += Math.min(pixels[0] + 75, 255);
                        pixels[1] += Math.min(pixels[1] + 10, 255);
                    }
                }



                raster.setPixel(xx, yy, pixels);
            }
        }
        return image;
    }
}
