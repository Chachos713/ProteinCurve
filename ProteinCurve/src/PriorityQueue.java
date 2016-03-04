/**
 * A simple implementation of a priority queue to help sort the atoms based on a
 * distance from a point.
 * 
 * @author Kyle Diller
 *
 * @param <E>
 *            the type to store in the queue.
 */
public class PriorityQueue<E> {
	/**
	 * The start of the queue of data to be stored.
	 */
	private Node<E> headData;

	/**
	 * The start of the queue of what the queue is sorted by.
	 */
	private Node<Double> headSort;

	/**
	 * The size of the queue.
	 */
	private int size;

	/**
	 * Creates a priority queue with nothing in it.
	 */
	public PriorityQueue() {
		headData = new Node<E>(null);
		headSort = new Node<Double>(null);
	}

	/**
	 * Adds a new data point to the queue. Uses a linear insertion method to
	 * determine where it fits.
	 * 
	 * @param data
	 *            the data to add to the queue.
	 * @param n
	 *            the value used to sort it by.
	 */
	public void add(E data, Double n) {
		size++;
		Node<E> tempData = headData;
		Node<E> prevData = headData;
		Node<Double> tempSort = headSort;
		Node<Double> prevSort = headSort;

		Node<E> d = new Node<E>(data);
		Node<Double> s = new Node<Double>(n);

		while (prevSort.hasNext()) {
			prevData = tempData;
			prevSort = tempSort;
			tempSort = tempSort.getNext();
			tempData = tempData.getNext();

			if (tempSort == null || tempSort.getData() > n) {
				break;
			}
		}

		d.setNext(prevData.getNext());
		prevData.setNext(d);
		s.setNext(prevSort.getNext());
		prevSort.setNext(s);
	}

	/**
	 * Check if there is still data in the queue to be read.
	 * 
	 * @return true if there is still data left in the queue.
	 */
	public boolean hasNext() {
		return headData.hasNext();
	}

	/**
	 * Gives the next piece of data in the queue.
	 * 
	 * @return The next piece of data in the queue.
	 */
	public E dequeue() {
		size--;
		Node<E> temp = headData.getNext();
		headSort.setNext(headSort.getNext().getNext());
		headData.setNext(temp.getNext());

		return temp.getData();
	}

	public String toString() {
		String n = "";

		Node<Double> d = headSort;
		Node<E> e = headData;

		while (d.hasNext()) {
			d = d.getNext();
			e = e.getNext();
			n += e.getData() + " <> " + d.getData() + "\n";
		}

		return n;
	}

	/**
	 * Gives the size of the queue.
	 * 
	 * @return the size of the queue.
	 */
	public int size() {
		return size;
	}
}
