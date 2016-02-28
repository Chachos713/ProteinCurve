public class Node<E> {
	private Node<E> next;
	private E data;

	public Node(E data) {
		this(data, null);
	}

	public Node(E data, Node<E> next) {
		this.data = data;
		this.next = next;
	}

	public E getData() {
		return data;
	}

	public void setData(E data) {
		this.data = data;
	}

	public Node<E> getNext() {
		return next;
	}

	public boolean hasNext() {
		return next != null;
	}

	public void setNext(Node<E> next) {
		this.next = next;
	}
}
