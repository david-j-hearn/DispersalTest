import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.PrintWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration;
//import java.ec.util.MersenneTwisterFast;

/**
 *
 * @author David
 */
public class calculateDispersalT {
		BufferedReader d = null;
        	String line;
        	String[] maxTaxa = new String[1000];
		String[] taxa = new String[1000];
        	boolean inTrees = false;
        	boolean inTaxa = false;
        	int cnt=0, cnt1=0;
        	boolean readingTrees = false;
		Hashtable hash;
		MersenneTwisterFast r = new MersenneTwisterFast();
		int[] numbers = null;

//Function options to add:
//	altnerate hypothesis: greater, lesser, or no directionality in time compared to null - RIGHT NOW IT IS SET AT LESSER
//	whether multiple morphological events are permissible at a node - RIGHT NOW, NO
//	format of the input file (newick, nexus, altnexus) - RIGHT NOW IT IS SET TO ALTNEXUS - BETTER TO SET TO NEWICK WITH A LIST OF THE TREES!!

public calculateDispersalT(String treeFile, String MorphologyMRCAFile, String DispersalMRCAFile, int nRandomizedReplicates, String outputFile)
	{
	try
		{
		String line1 = new String();
		PrintWriter out = new PrintWriter(new FileWriter(outputFile));
		//out.println("Tree number\tAverage distance between observed and randomized values for " + nRandomizedReplicates + "\tobserved average distance between events\tpercentation of replicates where observed is greater than random distance\tcount of replicates where observed in greater than random\n");
		out.println("Tree number\tRandomization replicate\tRandomized average distance between events\tobserved average distance between events\n");
			//out.println(nTrees + "\t" + (((double)totD)/nRandomizedReplicates) + "\t" + dispersalIndex + "\t" + Math.round(100*((double)cntSig/nRandomizedReplicates)) + "\t" + cntSig);

//read trees
		d = new BufferedReader(new FileReader(new File(treeFile)));

//foreach tree...
		double percentSig = 0;
		int nTrees=0;
		Tree tree = null;
		while((tree = this.treesFromAltnexusFileSingle())!=null)
		//while((tree = this.treesFromNexusFileSingle())!=null)
		//for(Tree tree: trees)
			{
			if(numbers==null)
				{
				taxa = Tree.getTaxonLabels(tree);
				numbers = new int[taxa.length-1];
				for(int i=0; i<numbers.length; i++) { numbers[i] = i; }
				}
//System.out.println("Top");
			nTrees++;
//System.out.println("On tree " + nTrees);
//System.out.println("Tree is " + Tree.treeToStringTimeOrder(tree,false));

//place the dispersal events
			int cntD = 0;
			BufferedReader br1 = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(DispersalMRCAFile))));
  			while ((line1 = br1.readLine()) != null)   
				{
				String[] data = line1.split("\\s+");
				try 
					{
					tree = Tree.setEventAtMRCA(tree, data[0],data[1], "dispersal");
					}
				catch(Exception e) 
					{
					System.out.println(data[0] + " and " + data[1] + " are the taxa whose MRCA we are finding.");
					System.out.println("Could not place dispersal event: " + e.getMessage());
					
					e.printStackTrace();
					System.exit(0);
					}
				cntD++;
  				}
  			br1.close();

			int cntM = 0;
//place the morphology transformation events
//System.out.println("M");
			double tot = 0.0;
			BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(MorphologyMRCAFile))));
  			while ((line1 = br.readLine()) != null)   
				{
				String[] data = line1.split("\\s+");
//System.out.println("The morphology placement taxa are " + data[0] + " and " + data[1]);
				try {
					tot += calculateDistanceToNearestOlderDispersalEvent(data[0],data[1],tree);
					}
				catch(Exception e) 
					{
					System.out.println("Could not place morphological event: " + e.getMessage());
					e.printStackTrace();
					System.exit(0);
					}
				cntM++;
  				}
  			br.close();

//calculate dispersal index
//System.out.println("D");
			//out.println("Observed dispersal distance:");
			double dispersalIndex = tot/cntM;
			//out.println(dispersalIndex);

//carry out randomizations
//System.out.println("R");
			//out.println("Randomized dispersal distances:");
			//double[] randomizedReps = new double[nReps];
			int cntSig = 0;
			double totD=0.0;
			for(int i=0; i<nRandomizedReplicates; i++)
				{
//if(i%1000 == 0) System.out.println("\t\tRandomization replicate " + i);

				//randomizedReps[i] = calculateRandomizedDistanceMetric(cntM, cntD, taxa);
				//double dist = calculateRandomizedDistanceMetric(cntM, cntD, tree);
				double dist = calculateRandomizedDistanceMetric(cntM, DispersalMRCAFile, tree);
				//out.println(nTrees + "\t" + i + "\t" + dist);
				//String yes = "";
				totD+= dist-dispersalIndex;
	out.println(nTrees + "\t" + i + "\t" + dist + "\t" + dispersalIndex);
				if(dist>dispersalIndex) 
					{
					cntSig++;
					//yes = "X";
					}
				//out.println(nTrees + "\t" + i + "\t" + dist + "\t" + dispersalIndex + "\t" + cntSig + "\t" + ((double)cntSig/(i+1)) + "\t" + yes);
			
				}
			totD /= (double)nRandomizedReplicates;
			System.out.println(nTrees + "\t" + totD);	
				
	//out.println(nTrees + "\t" + (((double)totD)/nRandomizedReplicates) + "\t" + dispersalIndex + "\t" + Math.round(100*((double)cntSig/nRandomizedReplicates)) + "\t" + cntSig);
//System.out.println(nTrees + "\t" + (((double)totD)/nRandomizedReplicates) + "\t" + dispersalIndex + "\t" + Math.round(100*((double)cntSig/nRandomizedReplicates)));

			if(Math.round(100*((double)cntSig/nRandomizedReplicates))>=95)
				{
				percentSig++;
				}
			
			if(nTrees%10 == 0) 
				{
				//System.out.println("Processing tree " + nTrees + ": " + percentSig + " of the trees have statistically significant times.");
				//System.out.println("Processing tree " + nTrees);
				out.flush();
				}
//System.out.println("Bottom");
			}
		d.close();
		if(nTrees ==0)
			{
			System.out.println("No trees were found in the file " + treeFile);
			System.exit(0);
			}
		percentSig /= nTrees;
		//out.println(percentSig + " of the " + nTrees + " trees had statistically significantly association between dispersal events and evolutionary transition events.");
		//System.out.println(percentSig + " of the " + nTrees + " trees had statistically significantly association between dispersal events and evolutionary transition events.");
		//System.out.println(cntSig + " of the " + nTrees + " trees had statistically significantly association between dispersal events and evolutionary transition events.");
		out.close();
		}
	catch(Exception p)
		{
		p.printStackTrace();
		}
	}


public int[] fisherYatesShuffle(int[] a)
	{
	int N = a.length;

        // shuffle
        for (int i = 0; i < N; i++) 
		{
            	// int from remainder of deck
            	int r = i + (int) (Math.random() * (N - i));
            	int swap = a[r];
            	a[r] = a[i];
            	a[i] = swap;
        	}
	return a;
	}


//public double calculateRandomizedDistanceMetric(int nMorphologyEvents, int nDispersalEvents, Tree tree)
public double calculateRandomizedDistanceMetric(int nMorphologyEvents, String DispersalMRCAFile, Tree tree)
	{
//clean tree
//System.out.println("In randomization");
//System.out.println("got taxa");
	tree = Tree.cleanTree(tree);
	double treeLen = Tree.getTotalLength(tree,true);
//System.out.println("The total tree length is " + treeLen);
//System.out.println(Tree.treeToStringTimeOrder(tree,false));

// PICK A NODE BETWEEN 0 AND N-1
//PLACE THE CONSTRAINT AT THE RANDOM NODE


//ALWAYS INCLUDE THE ROOT, AND ASSUME THAT THE USER PUT IN THE ROOT NODE AS WELL!!!
//ASSUMES PERFECTLY BIFURCATING TREE!!
/*
	tree = Tree.findRoot(tree);
	tree.event = "dispersal";

	for(int i=0; i<nDispersalEvents-1; i++)
		{
		double position = r.nextDouble()*treeLen;
		tree = Tree.findNode(tree,position);
		if(tree.event == null)
			{	
			tree.event = "dispersal";
			}
		else if(tree.event != null)
			{
			if(tree.event.equals("dispersal"))
				{
				i--;
				}
			}
		}
*/
			int cntD = 0;
			try{
			BufferedReader br1 = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(DispersalMRCAFile))));
			String line1;
  			while ((line1 = br1.readLine()) != null)   
				{
				String[] data = line1.split("\\s+");
				try 
					{
					tree = Tree.setEventAtMRCA(tree, data[0],data[1], "dispersal");
					}
				catch(Exception e) 
					{
					System.out.println(data[0] + " and " + data[1] + " are the taxa whose MRCA we are finding.");
					System.out.println("Could not place dispersal event: " + e.getMessage());
					
					e.printStackTrace();
					System.exit(0);
					}
				cntD++;
  				}
  			br1.close();
			}
			catch(Exception e1) 
				{
				e1.printStackTrace();
				System.exit(1);
				}

	double tot=0.0;
	double totAttempts = 0;
	for(int i=0; i<nMorphologyEvents; i++)
		{
		double position = r.nextDouble()*treeLen;
		tree = Tree.findNode(tree, position);
		if(!tree.checker)
			{
			tree = Tree.checkOffClade(tree);
			if(tree == null)
				{
				System.out.println("Placing morphology event failed.\nQuitting.\n");
				System.exit(0);
				}
			//tree.checker = true;		//morphology event already at this node!!
//trace backwards until dispersal event or root is found
			double time = 0.0;
			while(tree.ancestor!=null)
				{
				if(tree.event!=null)
				if(tree.event.equals("dispersal"))
					{
//System.out.println("Found dispersal event");
					break;
					}
				time+=tree.time;
//System.out.println("The node is " + tree.value + " The time is " + tree.time + " The total time is " + time);

				tree = tree.ancestor;
				}
			tot+=time;
//System.out.println("The node is " + tree.value + " and the event is " + tree.event + " and the total time is " + tot + " \n\n");
			}
		else
			{
			i--;
			totAttempts++;
			if(totAttempts>20)
				{
				tot = 0.0;
				tree = Tree.cleanChecker(tree);
				i=0;
				}
			}
		}
	
	return tot/nMorphologyEvents;
	}


//t1 and t2 are two taxa that define the MCRA where the event of interest occurred
//assumes the node times are set
public double calculateDistanceToNearestOlderDispersalEvent(String t1, String t2, Tree t) throws Exception
	{
        Tree tree = new Tree(t);
        Tree copy = new Tree(t);
//System.out.println(Tree.treeToStringTimeOrder(tree,false));
        Tree copy1 = Tree.findTaxon(copy, t1, true);
        tree = Tree.findTaxon(tree, t2, true);

//System.out.println("In calculateDistanceToNearestOlderDispersalEvent");
	if(copy1==null || tree==null) 
		{ 
		System.out.println("Could not find taxa"); 
		System.exit(0);
		}

	double time = 0.0;

        if(copy1.value.equals(tree.value))
                {
//System.out.println("Matched!");
		if(tree.discreteState != null)
		if(tree.discreteState.equals("morphology"))
			throw new Exception("Morphological event already placed at MRCA of " + t1 + " and " + t2); 
		tree.discreteState = "morphology";		
		while(tree.ancestor != null)
			{
//System.out.println(tree.value);
			if(tree.event.equals("dispersal"))
				{
//System.out.println(time);
				return time;
				}
			time += tree.time;
			tree = tree.ancestor;
			}
//System.out.println(time);
		return time;
                }

        while(tree.ancestor!=null)
                        {
                        copy = copy1;
//System.out.println(tree.value + ", " + copy.value);
                        tree = tree.ancestor;
//System.out.println(tree.value + ", " + copy.value);
                        while(copy.ancestor!=null)
                                {
                                copy = copy.ancestor;
//System.out.println(tree.value + ", " + copy.value);
                                if(copy.value.equals(tree.value))
                                        {
//System.out.println("Matched!");
					if(tree.discreteState != null)
					if(tree.discreteState.equals("morphology"))
						throw new Exception("Morphological event already placed at MRCA of " + t1 + " and " + t2); 
					tree.discreteState = "morphology";		
					while(tree.ancestor != null)
						{
//System.out.println(tree.value);
						if(tree.event.equals("dispersal"))
							{
//System.out.println(time);
							return time;
							}
						time += tree.time;
						tree = tree.ancestor;
						}
//System.out.println(time);
					return time;
                                        }
                                }
                        }
	throw new Exception("Could not find common ancestor for " + t1 + " and " + t2);
	}

//public Tree[] treesFromAltnexusFile()
public Tree treesFromAltnexusFileSingle()
        {
        try{

        String tree;
        String label = new String();
	//Vector<Tree> treeVector = new Vector();
        while((line=d.readLine())!=null)
                {
                if(inTrees && line.toLowerCase().matches(".*end;*")) 
			{
//System.out.println("A");
			//return (Tree[])treeVector.toArray(new Tree[treeVector.size()]);
			return null;
			}
                else if(inTrees && line.toLowerCase().matches(".*;.*"))
                        {
                        for(int k = 0; k<line.length(); k++)
                                {
                                char c = line.charAt(k);
                                if(c=='[')
                                        {
                                        while(line.charAt(k+1)!=']')
                                                {
                                                k++;
                                                }
                                        k++;
                                        }
                                else
                                        {
                                        label += line.charAt(k);
                                        }
                                }
                        String[] dta = label.split("\\s*=\\s*");
                        tree = dta[dta.length-1];
//System.out.println("Tree is " + tree);
			//treeVector.add(new Tree(tree,true));
                        return(new Tree(tree,true));
                        }
                if(line.toLowerCase().matches(".*begin trees;.*")) inTrees = true;
                else if(inTrees && line.toLowerCase().matches(".*end;*")) 
			{
//System.out.println("B");
			//return (Tree[])treeVector.toArray(new Tree[treeVector.size()]);
			return null;
			}
                }
        //d.close();
        }
        catch(Exception e){e.printStackTrace();}
	//return (Tree[])treeVector.toArray(new Tree[treeVector.size()]);
//System.out.println("C");
        return null;
        }

public Tree treesFromNexusFileSingle()
        {
        try{

        String tree;
        String label = new String();
	//Vector<Tree> treeVector = new Vector();
        while((line=d.readLine())!=null)
                {
                if(inTrees && inTaxa && !readingTrees)
                        {
                        if(line.toLowerCase().matches(".*;.*")) readingTrees = true;
                        line = line.replaceAll(",","");
                        line = line.replaceAll(";","");
                        line = line.replaceAll("^\\s+", "");
                        String[] data = line.split("\\s");
                        maxTaxa[Integer.parseInt(data[0])-1] = data[1];
                        cnt++;
                        }
                else if(readingTrees)
                        {
                        for(int k = 0; k<line.length(); k++)
                                {
                                char c = line.charAt(k);
                                if(c=='[')
                                        {
                                        while(line.charAt(k+1)!=']')
                                                {
                                                k++;
                                                }
                                        k++;
                                        }
                                else
                                        {
                                        label += line.charAt(k);
                                        }
                                }
                        String[] dta = label.split("\\s*=\\s*");
                        tree = dta[dta.length-1];
                        for(int j=0;j<cnt;j++)
                                {
                                tree = tree.replaceAll("\\(" + (j+1) + ":", "\\(" + maxTaxa[j] + ":");
                                tree = tree.replaceAll("\\(" + (j+1) + ",", "\\(" + maxTaxa[j] + ",");
                                tree = tree.replaceAll("," + (j+1) + ":", "," + maxTaxa[j] + ":");
                                tree = tree.replaceAll("," + (j+1) + ",", "," + maxTaxa[j] + ",");
                                tree = tree.replaceAll("," + (j+1) + "\\)", "," + maxTaxa[j] + "\\)");
                                }
			//treeVector.add(new Tree(tree,true));
                        return(new Tree(tree,true));
                        }
                if(line.toLowerCase().matches(".*begin trees;.*")) inTrees = true;
                else if(line.toLowerCase().matches(".*translate.*")) inTaxa = true;
                else if(inTrees && line.toLowerCase().matches(".*end;*")) 
			{
			//return (Tree[])treeVector.toArray(new Tree[treeVector.size()]);
			return null;
			}
                }
        //d.close();
        }
        catch(Exception e){e.printStackTrace();}
	//return (Tree[])treeVector.toArray(new Tree[treeVector.size()]);
        return null;
        }

}

