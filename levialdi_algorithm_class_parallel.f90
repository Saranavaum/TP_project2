program levialdi_algorithm

    use dll
    use mpi_f08
    
    implicit none
    
    !------------------------Variables--------------------------
    type(DoubleLinkedList):: imag_hist
    type(Node), pointer :: current_node
    integer :: m,n,n_rows,n_cols
    logical, dimension(:, :), pointer :: old_imag, new_imag ,tmp_imag, big_imag
    integer :: K=0
    integer, dimension(:,:), allocatable ::old_labels(:,:),new_labels(:,:)
    integer ::  row, col, ib, ie, jb, je   
    integer :: n_ranks, my_rank, root, north_rank, south_rank, east_rank, west_rank
    logical :: root_still
    type(MPI_Status) :: status
    type(MPI_Datatype) :: a_col,long_row,short_row
    
    
    
    
    
    !-------------------Program--------------------------
    !inicializamos mpi
    call MPI_Init()
    call MPI_Comm_size(MPI_COMM_WORLD, n_ranks)
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
    
    root=0
    !leemos los números de entrada en el proceso raíz (root rank) 
    if (my_rank == root) then
        read *, m, n, n_rows, n_cols
    
    end if
    ! We distribute to all processes using the MPI_Bcast function. We make sure that all processes have the necessary data to work.

    call MPI_Bcast(m, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(n, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(n_rows, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(n_cols, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    
    
    if (my_rank==root) then
        if (n_rows * n_cols /= n_ranks) then
            print "(a,i0,a)", "Run this code with `-np ", n_cols*n_rows , "`."
            call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_TOPOLOGY)
        end if
    end if
    
    !we define the coordinates of each rank    
    call get_coords( my_rank, n_rows, n_cols, row, col )
    north_rank = get_rank( row - 1, col,     n_rows, n_cols )
    south_rank = get_rank( row + 1, col,     n_rows, n_cols )
    west_rank  = get_rank( row,     col - 1, n_rows, n_cols )
    east_rank  = get_rank( row,     col + 1, n_rows, n_cols )
    
    !We calculate the indices of the partition of the matrix that will go to each rank in the most optimal way possible, whether it is even or odd.
    call partition(row, n_rows, m, ib, ie)
    call partition(col, n_cols, n, jb, je)
    
        
    allocate(old_imag(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(tmp_imag(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(new_imag(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(big_imag(m,n))

    allocate(old_labels(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(new_labels(ib - 1:ie + 1, jb - 1:je + 1))
    
        
    nullify(imag_hist%head)
    nullify(imag_hist%tail)
    
    
    !We define the MPI_Datatype for the rows and columns taking into account the indices of the partitions we did previously and the row stride.
    call MPI_Type_contiguous( ie - ib + 1, MPI_LOGICAL, a_col )
    call MPI_Type_commit( a_col )
    
    
    block
        type(MPI_Datatype) :: a_tmp_row, another_tmp_row
        integer(kind=MPI_ADDRESS_KIND) :: lb, real_extent

        call MPI_Type_vector( je - jb + 3, 1, ie - ib + 3, MPI_LOGICAL, a_tmp_row )
        call MPI_Type_get_extent( MPI_LOGICAL, lb, real_extent )
        call MPI_Type_create_resized( a_tmp_row, lb, real_extent, long_row )
        call MPI_Type_commit( long_row )
        
        call MPI_Type_vector( je - jb + 1, 1, ie - ib + 3, MPI_LOGICAL, another_tmp_row )
        call MPI_Type_get_extent( MPI_LOGICAL, lb, real_extent )
        call MPI_Type_create_resized( another_tmp_row, lb, real_extent, short_row)
        call MPI_Type_commit( short_row )
    end block
    
    !We read the image
    call read_map( old_imag, m, n )
    !We set the type of topology for the borders, which in this case will not be a ring, it's null
    call update_borders( old_imag)
    !We synchronize the borders with the rest of the ranks.
    call synchronize_ghost(old_imag)
    
    !We add to a list the fragment calculated by each rank.
    call add_node(imag_hist,old_imag)
    
    
    !We can join all the fragments and save the entire image in the following way, but it is a bit redundant and for the growing part it is more efficient to keep it fragmented:
    !call save_map(old_imag,m,n,big_imag)
    !call add_node(imag_hist,big_imag)
    
    
    !Shrinking components: 
    do 
        call update_pixels( old_imag, new_imag)
        call update_borders( new_imag)
        call synchronize_ghost(new_imag)
        
     !   call save_map(new_imag,m,n,big_imag)
        call add_node(imag_hist,new_imag)
        
        !all ranks communicate to know when all of them have reached a .false. value."
        root_still = world_is_still( new_imag)
        call MPI_Bcast(root_still, 1, MPI_LOGICAL, root, MPI_COMM_WORLD)
        if (root_still) exit

        tmp_imag => old_imag;  old_imag => new_imag;  new_imag => tmp_imag
    end do 
    

    !The code works correctly up to this point, but the part of the growing still needs to be implemented properly so that all ranks have updated information on the value of K, but I haven't been able to get there yet.
    
    
    !Checking Shrinking components works corretly:
   ! if (my_rank==1) then
    !    call printList(imag_hist%head )
     !   print*,'my rank is',my_rank
      !  print*,'old_antes'
    !end if     
    
    
    !growing components
    old_labels=0
    current_node => imag_hist%tail

    do while (associated(current_node))
      
        !call growing(current_node%myData,old_labels, new_labels,m,n,K)

       
       

       ! old_labels=new_labels
        if (associated(current_node, imag_hist%head)) exit
        current_node => current_node%prev
    end do
    
    
    
    
    

    call destroyList(imag_hist%head)
    if (associated( old_imag)) deallocate(old_imag)
    if (associated( new_imag)) deallocate(new_imag)

 
    call MPI_Type_free( a_col )
    call MPI_Type_free( short_row )
    call MPI_Type_free( long_row )
    call MPI_Finalize() 
contains


    function world_is_still( old_map) result(world_still)
        logical, dimension(:, :), pointer, intent(in) :: old_map
        logical,dimension(ib - 1:ie + 1, jb - 1:je + 1):: is_false
        logical :: world_still
        logical, dimension(n_ranks)                   :: all_is_still
        is_false=.False.
        world_still = all( old_map .eqv. is_false )
        call MPI_Gather(world_still, 1, MPI_LOGICAL, &
                        all_is_still,   1, MPI_LOGICAL, root, MPI_COMM_WORLD)
        if (my_rank == root) world_still = all ( all_is_still .eqv. .True.)

    end function world_is_still






    subroutine update_pixels( old_map, new_map)
        logical, dimension(ib - 1:ie + 1, jb - 1:je + 1), intent(inout) :: old_map
        logical, dimension(ib - 1:ie + 1, jb - 1:je + 1), intent(inout) :: new_map
        integer :: i, j
        
        do j = jb, je
            do i = ib, ie
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
      integer  ,dimension(ib - 1:ie + 1, jb - 1:je + 1),intent(inout):: old_label
      integer  ,dimension(ib - 1:ie + 1, jb - 1:je + 1),intent(inout):: new_label
      integer :: i, j, c,k
  
 
      do j = jb, je
          do i = ib, ie
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
    
    
    subroutine get_coords( rank, n_rows, n_cols, row, col )
        integer, intent(in)    :: rank, n_rows, n_cols
        integer, intent(inout) :: row, col

        row = modulo(rank, n_rows)
        col = (rank - row) / n_rows
        if (0 <= col .and. col < n_cols) then
            return
        else
            print "(a, 2(i0, a))", "get_coords: rank ", rank, &
                " is outside the column range [0, ", n_cols, ")."
            call MPI_Abort( MPI_COMM_WORLD, MPI_ERR_TOPOLOGY )
        end if
    end subroutine get_coords    
    
    
    
    subroutine partition( id, n_ids, size, b, e )
        integer, intent(in)    :: id, n_ids, size
        integer, intent(inout) :: b, e
        integer :: remainder, quotient

        remainder = modulo( size, n_ids )
        quotient  = (size - remainder) / n_ids
        b = 1 + quotient * (id    ) + min( remainder, id     )
        e =     quotient * (id + 1) + min( remainder, id + 1 )
    end subroutine partition    


    integer function get_rank( row, col, n_rows, n_cols )
        integer, intent(in) :: row, col, n_rows, n_cols

        if (      0 <= col .and. col < n_cols &
            .and. 0 <= row .and. row < n_rows) then
                get_rank = row + col * n_rows
        else
            get_rank = MPI_PROC_NULL
        end if
    end function get_rank


    subroutine read_map( map, h, w )
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: h, w
        character(len=:), allocatable :: line
        logical,          allocatable :: temp(:)
        integer :: i, j, rb, re, cb, ce
       
        if (my_rank == root) then
            block
                integer :: current_row
                integer :: current_col
                integer :: dst
                allocate(character(len=w) :: line)
                do current_row = 0, n_rows - 1
                    call partition(current_row, n_rows, h, rb, re)
                    do i = rb, re
                        read *, line
                        do current_col = 0, n_cols - 1
                            call partition(current_col, n_cols, w, cb, ce)
                            dst = get_rank(current_row, current_col, n_rows, n_cols)
                            allocate(temp(ce - cb + 1))
                            do j = cb, ce 
                                select case (line(j:j))
                                case ('X')
                                    temp(j - cb + 1) = .True.
                                case ('.')
                                    temp(j - cb + 1) = .False.
                                case default
                                    stop "read_map: wrong input character `" // line(j:j) // "`"
                                end select
                            end do
                            if (dst == root) then
                                map(i, cb : ce)  = temp
                            else
                                call MPI_Send( temp, ce - cb + 1, MPI_LOGICAL, dst, 0,  MPI_COMM_WORLD )
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                    end do
                end do
                if (allocated( line )) deallocate(line)
            end block
        else
            do i = ib, ie
                call MPI_Recv(map(i,jb), 1, short_row, root, 0, MPI_COMM_WORLD, status)
            end do
        end if
    end subroutine read_map




    subroutine save_map(map, h, w, big_map)
        logical, dimension(:, :), pointer, intent(in) :: map
        integer, intent(in)                           :: h, w
        logical, dimension(:, :), pointer, intent(inout) :: big_map
        
        logical,          allocatable :: temp(:)
        integer :: i, j, rb, re, cb, ce

        if (my_rank == root) then
            block
                integer :: current_row
                integer :: current_col
                integer :: dst

                do current_row = 0, n_rows - 1
                    call partition(current_row, n_rows, h, rb, re)
                    do i = rb, re
                        do current_col = 0, n_cols - 1
                            call partition(current_col, n_cols, w, cb, ce)
                            allocate(temp(ce - cb + 1))
                            dst = get_rank(current_row, current_col, n_rows, n_cols)
                            if (dst == root) then
                                do j = cb, ce 
                                    big_map(i, j) = map(i,j)
                                end do
                            else 
                                call MPI_Recv(temp, ce-cb+1, MPI_LOGICAL, dst, 0, MPI_COMM_WORLD, status)
                                do j = cb, ce 
                                    big_map(i, j) = temp(j-cb+1)
                                end do
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do

                        
                    end do
                end do
                print *
            end block
       else
           do i = ib, ie
               call MPI_Send( map(i, jb), 1, short_row, root, 0, MPI_COMM_WORLD)
           end do
       end if
        

    end subroutine save_map



    subroutine print_int(map, h, w)
        integer, dimension(:, :), pointer, intent(in) :: map
        integer, intent(in)                           :: h, w       
        integer,          allocatable :: temp(:)
        integer :: i, j, rb, re, cb, ce

        if (my_rank == root) then
            block
                integer :: current_row
                integer :: current_col
                integer :: dst
                integer,dimension(:), allocatable :: line
                allocate(line(1:w))

                do current_row = 0, n_rows - 1
                    call partition(current_row, n_rows, h, rb, re)
                    do i = rb, re
                        do current_col = 0, n_cols - 1
                            call partition(current_col, n_cols, w, cb, ce)
                            allocate(temp(ce - cb + 1))
                            dst = get_rank(current_row, current_col, n_rows, n_cols)
                            if (dst == root) then
                                do j = cb, ce 
                                    line(j:j) = map(i,j)
                                end do
                            else 
                                call MPI_Recv(temp, ce-cb+1, MPI_INTEGER, dst, 0, MPI_COMM_WORLD, status)
                                do j = cb, ce 
                                    line(j:j) = temp(j-cb+1)
                                end do
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                        print *, line
                        
                    end do
                end do
                print *
                if (allocated( line )) deallocate(line)
            end block
       else
           do i = ib, ie
               call MPI_Send( map(i, jb), 1, short_row, root, 0, MPI_COMM_WORLD)
           end do
       end if
        
    end subroutine print_int




    subroutine read_inter( map, h, w )
        logical, dimension(:, :),intent(inout) :: map
        integer, intent(in) :: h, w
        logical,          allocatable :: temp(:)
        integer :: i, rb, re, cb, ce
       
        if (my_rank == root) then
            block
                integer :: current_row
                integer :: current_col
                integer :: dst
                do current_row = 0, n_rows - 1
                    call partition(current_row, n_rows, h, rb, re)
                    do i = rb, re
                        do current_col = 0, n_cols - 1
                            call partition(current_col, n_cols, w, cb, ce)
                            dst = get_rank(current_row, current_col, n_rows, n_cols)
                            allocate(temp(ce - cb + 1))
                            if (dst /= root) then
                                call MPI_Send( temp, ce - cb + 1, MPI_integer, dst, 0,  MPI_COMM_WORLD )
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                    end do
                end do
            end block
        else
            do i = ib, ie
                call MPI_Recv(map(i,jb), 1, long_row, root, 0, MPI_COMM_WORLD, status)
            end do
        end if
    end subroutine read_inter





    subroutine synchronize_ghost( map )
        logical, dimension(:, :), pointer, intent(inout) :: map

        ! syncrhonize cols 
        call MPI_Sendrecv( map(ib, jb), 1, a_col, west_rank, 1, &
                           map(ib, je+1), 1, a_col, east_rank, 1, MPI_COMM_WORLD, status)
        call MPI_Sendrecv( map(ib, je), 1, a_col, east_rank, 2, &
                           map(ib, jb - 1), 1, a_col, west_rank, 2, MPI_COMM_WORLD, status)
        

        ! synchronize rows (big rows)
        call MPI_Sendrecv( map(ib,jb-1),     1, long_row, north_rank, 3, &
                           map(ie + 1,jb-1), 1, long_row, south_rank, 3, MPI_COMM_WORLD, status)
        call MPI_Sendrecv( map(ie,jb-1),     1, long_row, south_rank, 4, &
                           map(ib-1, jb-1),   1, long_row, north_rank, 4, MPI_COMM_WORLD, status)


    end subroutine synchronize_ghost



    subroutine update_borders( old_mtx)
        logical, dimension(:, :), pointer, intent(in out)     :: old_mtx

        ! Inner columns:
        old_mtx(ib:ie, jb - 1) = .False.
        old_mtx(ib:ie, je + 1) = .False.
        ! Full rows:
        old_mtx(ib - 1, jb - 1:je + 1) = .False.
        old_mtx(ie + 1, jb - 1:je + 1) = .False.
    end subroutine update_borders



end program levialdi_algorithm
