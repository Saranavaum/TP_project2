module dll
  implicit none

  type :: Node
     logical, dimension(:,:),allocatable  :: myData !cambiar esto como un puntero
     type(Node), pointer :: next
     type(Node), pointer :: prev !esto no hace falta
  end type Node

  type :: DoubleLinkedList
     type(Node), pointer :: head
     type(Node), pointer :: tail
  end type DoubleLinkedList

  contains

  subroutine add_node(list, matrix)
    type(DoubleLinkedList), intent(inout) :: list
    logical, dimension(:,:), intent(in) :: matrix

    integer :: height
    integer :: width
    ! Create new node
    type(Node), pointer :: new_node

    height = size(matrix, 1)
    width = size(matrix, 2)

    allocate(new_node)
    allocate(new_node%myData(height , width  ))



    new_node%myData(1: height, 1: width) = matrix(:,:) !esto lleva mucho tiempo,es mejor que solo copie de memory addres
    new_node%next => null()
    new_node%prev => list%tail

    ! Add new node to list
    if (associated(list%tail)) then
       list%tail%next => new_node
    else
       list%head => new_node
    end if
    list%tail => new_node

  end subroutine add_node
 
  recursive subroutine printList(head)
    type(Node), pointer, intent(in) :: head
    
    if(associated(head)) then
      call print_matrix(head%myData)
      call printList(head%next)
    end if 
  end subroutine printList

  subroutine print_matrix(matrix)
    logical, dimension(:,:), intent(in) :: matrix
    integer :: i, j

    
    do i = 1, size(matrix, 1)
      do j = 1, size(matrix, 2)
        if (matrix(i,j)) then
          write (*, "(A1)", advance="no") "X "
        else
          write (*, "(A1)", advance="no") ". "
        end if
      end do
      print *
    end do
    print *
    
  end subroutine print_matrix
 
  recursive subroutine destroyList(head)
    type(node), pointer, intent(inout) :: head
    
    if (associated(head)) then
      call destroyList(head%next)
      deallocate(head%myData)
      deallocate(head)
    end if
    
  end subroutine destroyList

end module dll 
