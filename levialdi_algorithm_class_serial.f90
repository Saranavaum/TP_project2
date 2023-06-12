program levialdi_algorithm
    use dll
    implicit none
    type(DoubleLinkedList):: imag_hist
    type(Node), pointer :: current_node
    integer :: m,n,nrows,ncols
    logical, dimension(:, :), pointer :: old_imag, new_imag ,tmp_imag
    integer :: K=0,i,j
    integer, dimension(:,:), pointer::old_labels(:,:),new_labels(:,:)
    
    read *, m, n, nrows,ncols 
    allocate(old_imag(m,n))
    allocate(tmp_imag(m,n))
    allocate(new_imag(m,n))
    

    allocate(old_labels(0:m+1,0:n+1))
    allocate(new_labels(0:m+1,0:n+1))
    
    nullify(imag_hist%head)
    nullify(imag_hist%tail)
    
    call read_map( old_imag, m, n )
    
    call add_node(imag_hist,old_imag)
    

    
    !Shrinking component
    do while(.not.check_all_false(imag_hist%tail%myData))
        call update_pixels( imag_hist%tail%myData, new_imag, m, n )
        call add_node(imag_hist,new_imag)

        !tmp_imag => old_imag;  old_imag => new_imag;  new_imag => tmp_imag
    end do 
    
    !call printList(imag_hist%head )
    
    old_labels=0
    !growing componentes
    current_node => imag_hist%tail
    do while (associated(current_node))
        !call print_matrix(current_node%myData)
        
        call growing(current_node%myData,old_labels, new_labels,m,n,K)
        old_labels=new_labels
        
        if (associated(current_node, imag_hist%head)) exit
        current_node => current_node%prev
    end do
    
    
    print"(I0)",K
    do i = 1,m 
        do j = 1, n
            write(*, "(I4)", advance='no') new_labels(i,j)
        end do
        write(*,*)
    end do

    print"(I0)",K 


    
    call destroyList(imag_hist%head)
    if (associated( old_imag)) deallocate(old_imag)
    if (associated( new_imag)) deallocate(new_imag)
    if (associated( tmp_imag)) deallocate(tmp_imag)
    if (associated( new_labels)) deallocate(new_labels)
    if (associated( old_labels)) deallocate(old_labels)

contains




    subroutine read_map( map, h, w )
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: h, w
        character(len=:), allocatable :: line
        integer :: i, j
        
        allocate(character(len=w) :: line)
        do i = 1, h
            read *, line
            do j = 1, w
                select case (line(j:j))
                case ('X')
                    map(i, j) = .true.
                case ('.')
                    map(i, j) = .false.
                case default
                    stop "read_map: wrong input character `" // line(j:j) // "`"
                end select
            end do
        end do
        if (allocated( line )) deallocate(line)
    end subroutine read_map





    subroutine update_pixels( old_map, new_map, h, w )
        integer, intent(in) :: h, w
        logical, dimension(0:h+1, 0:w+1), intent(inout) :: old_map
        logical, dimension(:, :), intent(inout) :: new_map
        integer :: i, j
	

        do j = 1, w
            do i = 1, h
                !if (j==1 .and. i==2) then
               !     print*, old_map(i,j),old_map(i-1,j),old_map(i,j-1),old_map(i-1,j-1)
               ! end if
                if (old_map(i, j)) then ! is white
                    if (.not. old_map(i-1,j) .and. .not. old_map(i,j-1) .and. .not. old_map(i-1,j-1)) then
                        new_map(i,j)=.false.
                     else
                        new_map(i,j)=old_map(i,j)   
                    end if
                else ! cell is black
                    if ( old_map(i-1,j) .and. old_map(i,j-1).and. old_map(i-1,j-1)) then
                        new_map(i,j)=.true.
                     else if ( old_map(i-1,j) .and. old_map(i,j-1) .and. .not. old_map(i-1,j-1)) then
                        new_map(i,j)=.true.
                     else
                        new_map(i,j)=old_map(i,j)  
                    end if
                end if
                
            end do
        end do
    end subroutine update_pixels
    

        
    subroutine growing(h_map,old_label, new_label,h,w,k)
      integer, intent(in) :: h,w
      logical, dimension(0:h+1, 0:w+1),intent(inout) :: h_map
      integer  ,dimension(0:h+1, 0:w+1),intent(inout):: old_label
      integer  ,dimension(0:h+1, 0:w+1),intent(inout):: new_label
      integer :: i, j, c,k
  
      
      
      
      do j = 1, w
          do i = 1, h
              c = count( h_map(i - 1:i + 1, j - 1:j + 1))
              if (h_map(i,j)) then
                  !labels(i,j)=1
                  if (c==1) then
                      k=k+1
                      new_label(i,j)=k
                      
                  else
                      if (h_map(i+1,j) .and. old_label(i+1,j)/=0) then
                          new_label(i,j)=old_label(i+1,j)
                          continue
                      else if (h_map(i,j+1).and. old_label(i,j+1)/=0) then
                          new_label(i,j)=old_label(i,j+1)
                          continue
                      else if (h_map(i+1,j+1) .and. old_label(i+1,j+1)/=0) then
                          new_label(i,j)=old_label(i+1,j+1)
                          continue
                      else if (h_map(i-1,j-1).and. old_label(i-1,j-1)/=0) then
                          new_label(i,j)=old_label(i-1,j-1)
                          continue
                      else if (h_map(i-1,j) .and. old_label(i-1,j)/=0) then
                          new_label(i,j)=old_label(i-1,j)
                          continue
                      else if (h_map(i,j-1) .and. old_label(i,j-1)/=0) then
                          new_label(i,j)=old_label(i,j-1)
                          continue
                      else if (h_map(i-1,j+1).and. old_label(i-1,j+1)/=0) then
                          new_label(i,j)=old_label(i-1,j+1)
                          continue
                      else if (h_map(i+1,j-1) .and. old_label(i+1,j-1)/=0) then
                          new_label(i,j)=old_label(i+1,j-1)
                          continue
                      else if (h_map(i,j) .and. old_label(i,j)/=0) then
                          new_label(i,j)=old_label(i,j)
                      end if
                       
  
                      
                      
                  end if
              else       
                  new_label(i,j)=0  
              end if
          end do
      end do
      

    end subroutine growing
    
    
    function check_all_false(map_check) result(check)
       logical,dimension(:,:),intent(inout)::map_check
       logical::check
       integer::c
    	
       c=count(map_check)
       if (c>0) then
           check=.false.
       else
           check=.true.
       end if


    end function check_all_false

end program levialdi_algorithm
