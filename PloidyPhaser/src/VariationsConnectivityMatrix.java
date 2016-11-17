import java.util.ArrayList;
import java.util.List;



public class VariationsConnectivityMatrix {
	ReadList[][]varExpMat;
	int matSize;
	
	public VariationsConnectivityMatrix (int d) {
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
	

	public void printMatrix(){
		
		for (int r=0;r<matSize;r++){
			for (int c=0;c<matSize;c++){
				getReadList( r, c).printReadList();
			}
		}
	}
	
	public void printMatrix(int limit){
		for (int r=0;r<limit;r++){
			for (int c=0;c<limit;c++){
				System.out.println (" ReadList"+r+"/"+c+":");
				getReadList( r, c).printReadList();
			}
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
