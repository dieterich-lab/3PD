/**
 * 
 */
package primerDesign.dsc.nio;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.nio.IntBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

import primerDesign.util.MyExtendedMath;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * Encapsulates a MemoryMappedFile (~'file-based array') to store integers.
 * 
 * @author Sebastian Fršhler
 *
 */
public class MemoryMappedIntFile implements Serializable{
	private static final long serialVersionUID = 1L;
	private static NumberFormat format = NumberFormat.getInstance();
	private IntBuffer buffer;
	private int capacity;
	private int size;
	private File filename;
	private MappedByteBuffer mapping;
	private static final int MAX_CONTIG_SIZE = MyExtendedMath.floor(Integer.MAX_VALUE*0.74/4);  // the empirically determined max size of a memory mapped int file (by me!, on my mac pro)
	
	/**
	 * Initializes a new MemoryMappedIntFile of size 'size'.
	 * 
	 * @param filename the file to map in memory
	 * @param size the size of the file to map in memory (the number of integer values to map)
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public MemoryMappedIntFile(File filename, int size) throws FileNotFoundException, IOException{
		if(size > MAX_CONTIG_SIZE) throw new IllegalArgumentException("This datastructure can handle at most " + MAX_CONTIG_SIZE + " entries!");
		this.mapping = new RandomAccessFile(filename, "rw").getChannel().map (FileChannel.MapMode.READ_WRITE, 0, size*4);
		this.buffer = this.mapping.asIntBuffer();
		//this.buffer = new RandomAccessFile(filename, "rw").getChannel().map (FileChannel.MapMode.READ_WRITE, 0, size*4).asIntBuffer();
		this.capacity = size;
		this.size = size * 4;
		this.filename = filename;
	}
	
	/**
	 * Opens a pre-existing MemoryMappedIntFile.
	 * 
	 * @param filename the file to open and map in memory
	 * 
	 * @throws IOException
	 */
	public MemoryMappedIntFile(File filename) throws IOException{
        FileChannel channel = new RandomAccessFile (filename, "rw").getChannel();
        this.buffer = channel.map (FileChannel.MapMode.READ_WRITE, 0, channel.size()).asIntBuffer();
        this.capacity = (int) channel.size()/4;
        this.size = (int) channel.size();
        this.filename = filename;
	}
	
	public void put(int index, int value){
		if(index >= this.capacity) throw new IllegalArgumentException("The max index of this memory mapped int file is " + format.format(this.capacity));
		this.buffer.put(index, value);
	}

	public void putQuick(int index, int value){
		this.buffer.put(index, value);
	}
	
	public int get(int index){
		if(index >= this.capacity) throw new IllegalArgumentException("The max index of this memory mapped int file is " + format.format(this.capacity));
		return this.buffer.get(index);
	}
	
	public int getQuick(int index){
		return this.buffer.get(index);
	}
	
	public int length(){
		return this.capacity;
	}
	
	private void writeObject(java.io.ObjectOutputStream out) throws IOException{
		out.writeObject(new MMIntFileState(this.capacity, this.size, this.filename));
	}
	
	public boolean isInMemory(){
		return this.mapping.isLoaded();
	}
	
	private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException{
		MMIntFileState state = (MMIntFileState) in.readObject();
		this.filename = state.getFile();
		this.capacity = state.getCapacity();
		this.size = state.getSize();
		FileChannel channel = new RandomAccessFile (filename, "rw").getChannel();
        this.buffer = channel.map (FileChannel.MapMode.READ_WRITE, 0, channel.size()).asIntBuffer();
	}
	
	public static int getMaxContigSize(){
		return MemoryMappedIntFile.MAX_CONTIG_SIZE;
	}
	
	public void sort(int from, int to, Comparator comparator){
		Integer[] tmp = new Integer[to + 1];
		for(int i=0; i<=to; i++){
			tmp[i] = this.buffer.get(i);
		}
		Arrays.sort(tmp, from, to+1, comparator);
		for(int i=0; i<tmp.length; i++){
			this.buffer.put(i, tmp[i]);
		}
	}
	
	public static void main(String[] args) throws IOException, IOException{    
        SimpleTimer timer = new SimpleTimer();
        Random random = new Random();
        MemoryMappedIntFile mappedFile;
        ObjectArrayList indices = new ObjectArrayList();
        NumberFormat format = NumberFormat.getInstance();
        
        int subIndices = Integer.parseInt(args[0]);
        int queries = Integer.parseInt(args[1]);
        
        int index;
        int position;
        
        long meanInit = 0;
        long queryTime;
        
        String filenameStub = "/tmp/SFtempfile";
        int size = MemoryMappedIntFile.getMaxContigSize();
        
        for(int i=0; i<subIndices; i++){
        	//System.out.print("Init memory mapped file: " + i + " of size: " + size);
        	indices.add(mappedFile = new MemoryMappedIntFile(new File(filenameStub + "_" + i + ".tmp"), size));
        	//System.out.println(" - done in " + timer.getTimeString());
        	//System.out.println("Is this buffer completely held in RAM?: " + mappedFile.isInMemory());
        }
        
        for(int i=0; i<subIndices; i++){
        	//System.out.print("Init sub-index: " + i);
        	mappedFile = (MemoryMappedIntFile)indices.get(i);
        	for(int j=0; j<size; j++){
        		mappedFile.put(j, j);
        	}
        	long time = timer.getTime();
        	meanInit += time;
        	//System.out.println(" - done in " + format.format(time) + " ms");
        }
        
        
        //System.out.print("Randomly retrieving entries from either of the memory mapped files");
        for(int i=0; i<queries; i++){
        	index = random.nextInt(subIndices);
        	position = random.nextInt(size);
        	((MemoryMappedIntFile)indices.get(index)).get(position);
        }
        int time = (int) timer.getTime();
        queryTime = time;
        //System.out.println(" - done in " + format.format(time) + " ms");
        
        
        //System.out.println("---");
        System.out.println("SubIndices: " + subIndices + " queries: " + queries);
        System.out.println("Mean init time: " + format.format(((double)meanInit / subIndices)) + " ms");
        System.out.println("Avg queries per second: " + format.format(((double)queries / ((double)queryTime / 1000))));
        System.out.println("Mean query: " + queryTime);
        System.out.println("---");
		
	}
	
	private class MMIntFileState implements Serializable{
		int capacity;
		int size;
		File filename;
		
		public MMIntFileState(int capacity, int size, File filename){
			this.capacity = capacity;
			this.size = size;
			this.filename = filename;
		}
		
		public int getCapacity(){
			return this.capacity;
		}
		
		public int getSize(){
			return this.size;
		}
		
		public File getFile(){
			return this.filename;
		}
	}
}
