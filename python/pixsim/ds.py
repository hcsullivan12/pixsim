import numpy as np

class Node(object):
    """Simple node structure."""
    def __init__(self, 
                 data=None, 
                 next_node=None):
        self.data = data
        self.next_node = next_node

class NodeBST(object):
    """Simple node structure for BST."""
    def __init__(self, 
                 key=None,
                 data=None):
        self.key, self.data = key, data
        self.left, self.right  = None, None
        self.parent = None

    def update(self, data):
        np.append(self.data, data, axis=0)
        return

class LinkedList(object):
    """Simple linked list structure."""
    def __init__(self, head=None):
        self.head = head
        self.size = 0
        if self.head is not None:
            self.size = 1

    def insert(self, data):
        node = Node(data=data, next_node=self.head)
        self.head = node
        self.size += 1

    def __len__(self):
        return self.size

    def array(self):
        # this list is 'backwards', make it forwards
        nel = self.__len__()
        ret = np.zeros((nel,self.head.data.shape[0]))
        node = self.head
        count = 0
        while node is not None:
            ret[nel-1-count] = node.data
            node = node.next_node
            count += 1
        return ret

class BST(object):
    def __init__(self):
        self.root = None
        self.size = 0
        return

    def __len__(self):
        return self.size

    def insert(self, key, data):
        node = NodeBST(key=key, data=data)
        prt = None
        x = self.root
        while x is not None:
            prt = x 
            if node.key == x.key:
                x.update(data)
                return
            if node.key < x.key:
                x = x.left
            else:
                x = x.right
        node.parent = prt
        self.size += 1
        if prt is None:
            self.root = node
        elif key < node.parent.key:
            node.parent.left = node
        else:
            node.parent.right = node


