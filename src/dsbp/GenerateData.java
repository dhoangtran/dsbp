package dsbp;
import java.io.IOException;
import java.util.Random;

import dsbp.model.*;
import dsbp.util.*;

@SuppressWarnings("unused")
public class GenerateData {
	public static void main(String[] args) throws IOException {
		for (int i = 1; i <= 10; i++) {
			DataGenerator gen = new DataGenerator(4, 6, 6, 500 * i, 500 * i);
			gen.generate();
			String filePath = String.format("_../../data/in/36x4_%d.txt", i);
			gen.write(filePath);
		}
		
	}		
}
