/*
   ] * To change this template, choose Tools | Templates
 6 and open the template in the editor.
 */


import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.Random;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration;

/**
 *
 * @author David
 */
public class getDispersalIndex  {
public static void main(String[] args)
	{
	if(args.length != 5)
		{
		System.out.println("\n\nThis program calculates the average distance between when a morphological transition occurs and when a dispersal event occurs that precedes the morphological transition.\nThe program will also estimate randomized values\nUsage:\n\njava getDispersalIndex <file with altnexus-format tree> <file that specifies locations of morphological transition in tree> <file that specifies locations of dispersal events in tree> <number of randomizations> <output file>\n\nThe format of the files is as follows:\n\ntaxon1 taxon2\ntaxon3 taxon4\n...\n\nwhere the pairs of taxa on each line define an ancestral node - their most recent common ancestor - where the event (morphological transition or dispersal) of interest occurred.\n\nALWAYS INCLUDE THE ROOT AS A DISPERSAL EVENT!!\n\n");
		System.exit(0);
		}
	String treeFile = args[0];
	String dispersalMRCAFile = args[2];
	String morphologyMRCAFile = args[1];
	int nRandomizationReplicates = Integer.parseInt(args[3]);
	String outputFile = args[4];
	
	calculateDispersalT cdt = new calculateDispersalT(treeFile, morphologyMRCAFile, dispersalMRCAFile, nRandomizationReplicates, outputFile);
	}
}




