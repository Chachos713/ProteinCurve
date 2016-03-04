/**
 * A simple implementation of a node, used in the priority queue.
 * 
 * @author Kyle Diller
 *
 * @param <E>
 *            Any value that is used by the priority queue.
 */
public class Node<E> {
	/**
	 * The next node in the que.
	 */
	private Node<E> next;

	/**
	 * The data stored in the que.
	 */
	private E data;

	/**
	 * Creates a node with data, but no node after it.
	 * 
	 * @param data
	 *            the data to be stored in the node.
	 */
	public Node(E data) {
		this(data, null);
	}

	/**
	 * Creates a node that has a node after it and has data.
	 * 
	 * @param data
	 *            the data to be stored in the node.
	 * @param next
	 *            the node that follows this one.
	 */
	public Node(E data, Node<E> next) {
		this.data = data;
		this.next = next;
	}

	/**
	 * Gets the data in the node.
	 * 
	 * @return the data stored within this node.
	 */
	public E getData() {
		return data;
	}

	/**
	 * Changes the data stored within the node.
	 * 
	 * @param data
	 *            the new data to be stored in the node.
	 */
	public void setData(E data) {
		this.data = data;
	}

	/**
	 * Get the node that follow this one.
	 * 
	 * @return the next node in the list.
	 */
	public Node<E> getNext() {
		return next;
	}

	/**
	 * Checks if there is a node after it.
	 * 
	 * @return true if there is a node that follows this one.
	 */
	public boolean hasNext() {
		return next != null;
	}

	/**
	 * Changes who the next node is.
	 * 
	 * @param next
	 *            the new node to follow this one.
	 */
	public void setNext(Node<E> next) {
		this.next = next;
	}
}
