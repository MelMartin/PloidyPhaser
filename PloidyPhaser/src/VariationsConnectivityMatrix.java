import java.util.ArrayList;
import java.util.List;



public class VariationsConnectivityMatrix {
	ReadList[][]varExpMat;
	int matSize;
	VariationsManager varMan;
	
	public VariationsConnectivityMatrix (VariationsManager vm,int d) {
		varMan=vm;
		matSize=d;
		varExpMat=new ReadList[d][d];
		for (int r=0;r<d;r++){
			for (int c=0;c<d;c++){
				varExpMat[r][c]=new ReadList();
			}
		}
	}
	
	public ReadList getReadList(int r,int c){
		//System.out.println (" ReadList: "+r+"/"+c);
		return varExpMat[r][c];
	}

	public void addConection(int row, int col, Integer readIndex) {
		varExpMat[row][col].add(readIndex);
	}
	

	public void printMatrix(int limit){
		for (int r=0;r<limit;r++){
			for (int c=0;c<limit;c++){
				System.out.println (" ReadList"+r+"/"+c+":");
				getReadList( r, c).printReadList();
			}
		}
	}
	
	
	public void findSeedNodes(int minReadsThreshold) {
		int conectivityWeight;//number of Reads that connect that node with all its direct neighbours
		int connectivityDegree;//number of direct neighbour nodes connected to this node
		int prevVarConWeight=0;//keeps track of the number of reads (coneectivity weight) in the previous Variation
		VariationData vd;//current Variation (Node) being examined
		VariationData nextVd=varMan.varExprIds.get(varMan.varExpMatIndexes.get(0));//keeps track of the next Variation to the current one
		int totalWeightPerPos=0;//sum of nb or reads per Variation expressions at the same position(2 or 3 potential expressions depending on the Variation).		
		int totalConnectivityPerPos=0;//sum of connected nodes 5connectivity degree)  per Variation expressions  (at same position)		
		int nbOfVarsPerPos=0;//nb of different potential expressions per Position (2 or 3  depending on the Variation).
		int contigPloidy=(int) PloidyPhaser.ploidies.get(nextVd.name);//get the contig ploidy from first vcf call
		boolean changesPos=false;//if the following variation changes position
		
		//variables for dealing with indels positions
		/*
		boolean inDelPrecedes=false;
		int nbOfInDelsPosRemaining=0;
		int preIndelPos=0;//keeps track of the Variation Position that contains the indel
		int estimatedPloidyRemaining=contigPloidy;
		*/
		
		//count number of conexions per variation
		for (int r=0;r<matSize;r++){//for all variations (per row)
			conectivityWeight=0;
			connectivityDegree=0;
			int c;//column number
			vd=nextVd;
			if(r<(matSize-1))nextVd=varMan.varExprIds.get(varMan.varExpMatIndexes.get(r+1));
			System.out.print(" "+vd.pos+"/"+vd.expSignature+": ");

			for ( c=0;c<r;c++){//with the fixed row, move through the columns until r=c
				if(getReadList(r, c).size()>minReadsThreshold) {
					conectivityWeight+=getReadList(r, c).size();
					connectivityDegree++;
					//System.out.print(" "+r+"/"+c+":"+vcm.getReadList(r, c).size());
				}
			}
			c=r;//now we fix the column and move through the rest of the rows 
			for (int rr=c;rr<matSize;rr++){
				if(getReadList(rr, c).size()>minReadsThreshold) {
					conectivityWeight+=getReadList(rr, c).size();
					connectivityDegree++;
				}
			}
			System.out.print("w:"+conectivityWeight+" d:"+connectivityDegree);
			
			//sum of weights and connectivity degrees per position
			if (!changesPos){//if the next variation has the same position than this one 
				totalWeightPerPos+=conectivityWeight;
				totalConnectivityPerPos+=connectivityDegree;
			}else {//it is the first expression of the variation (REference)
				totalWeightPerPos=conectivityWeight;
				totalConnectivityPerPos=connectivityDegree;
				changesPos=false;
				//prevVarConWeight=conectivityWeight;
			}
			nbOfVarsPerPos++;
	
			//variables to estimate if a node is crossed by an isolated colour or it's a mix of diffferent haplotypes
			int expectedIsolatedNodeRatio=0;;//Theoretical ratio of totalNbOfReads/k , which correspond to an isolated colour in the connectivity graph
			int upLimitVariation=0;
			int downLimitVariaiton=0;
			
			
			

			if (nextVd.pos!=vd.pos  || r==(matSize-1)){
				changesPos=true;//if the next variation changes position, actualize changePos for next loop				
				System.out.println(" totW:"+totalWeightPerPos+" totD:"+totalConnectivityPerPos+" nbOfVarsPerPos:"+nbOfVarsPerPos);
				nbOfVarsPerPos=0;
			}else System.out.println();
			
			
			
		}
		
		
		
		
		
	}
	
	
	
    public class ReadList{ 
    	
    	List <Integer> ReadsIndexes;
    	
    	private ReadList (){
    		ReadsIndexes=new ArrayList<Integer>();
    	}
    	
    	public int size(){
    		return ReadsIndexes.size();
    	}
    	
    	private boolean isNull(){    		
    		return ReadsIndexes==null;
    	}
    	
    	private void initialize(){
    		ReadsIndexes=new ArrayList<Integer>();
    	}
    	
    	private boolean add(Integer readIndex){
    		return ReadsIndexes.add(readIndex);
    	}
    	
    	private void printReadList(){
    		for (int i=0;i<ReadsIndexes.size();i++){
				System.out.print (" "+ReadsIndexes.get(i));
			}
    		System.out.println ();
    	}
    	
    }
	
}
