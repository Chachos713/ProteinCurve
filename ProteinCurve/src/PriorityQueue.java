public class PriorityQueue<E> {
	private Node<E> headData;

	private Node<Double> headSort;

	private int size;

	public PriorityQueue() {
		headData = new Node<E>(null);
		headSort = new Node<Double>(null);
	}

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

	public boolean hasNext() {
		return headData.hasNext();
	}

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

	public int size() {
		return size;
	}
}
