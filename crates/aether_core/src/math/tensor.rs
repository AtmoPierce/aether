// tensors.rs - minimal tensor data structures (no math, no allocation)
#![allow(non_camel_case_types)]


use core::ops::{Index, IndexMut};

// -----------------------------------------------------------------------------
// Fixed-size rank-3 tensor: Tensor3<T, D0, D1, D2>
// -----------------------------------------------------------------------------
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tensor3<T, const D0: usize, const D1: usize, const D2: usize> { pub data: [[[T; D2]; D1]; D0] }
impl<T: Copy, const D0: usize, const D1: usize, const D2: usize> Tensor3<T, D0, D1, D2> {
pub const fn new(data: [[[T; D2]; D1]; D0]) -> Self { Self { data } }
pub const fn dims() -> (usize, usize, usize) { (D0, D1, D2) }
pub fn as_flat(&self) -> &[T] { unsafe { core::slice::from_raw_parts(self.data.as_ptr() as *const T, D0*D1*D2) } }
pub fn as_flat_mut(&mut self) -> &mut [T] { unsafe { core::slice::from_raw_parts_mut(self.data.as_mut_ptr() as *mut T, D0*D1*D2) } }
}
impl<T, const D0: usize, const D1: usize, const D2: usize> Index<(usize,usize,usize)> for Tensor3<T, D0, D1, D2> { type Output = T; fn index(&self, ijk:(usize,usize,usize))->&T{let(i,j,k)=ijk;&self.data[i][j][k]} }
impl<T, const D0: usize, const D1: usize, const D2: usize> IndexMut<(usize,usize,usize)> for Tensor3<T, D0, D1, D2> { fn index_mut(&mut self, ijk:(usize,usize,usize))->&mut T{let(i,j,k)=ijk;&mut self.data[i][j][k]} }


// -----------------------------------------------------------------------------
// Fixed-size rank-4 tensor: Tensor4<T, D0, D1, D2, D3>
// -----------------------------------------------------------------------------
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tensor4<T, const D0: usize, const D1: usize, const D2: usize, const D3: usize> { pub data: [[[[T; D3]; D2]; D1]; D0] }
    impl<T: Copy, const D0: usize, const D1: usize, const D2: usize, const D3: usize> Tensor4<T, D0, D1, D2, D3> {
        pub const fn new(data: [[[[T; D3]; D2]; D1]; D0]) -> Self { 
            Self { 
                data 
            } 
        }
        pub const fn dims() -> (usize, usize, usize, usize) { 
            (D0, D1, D2, D3) 
        }
    }
    impl<T, const D0: usize, const D1: usize, const D2: usize, const D3: usize> Index<(usize,usize,usize,usize)> for Tensor4<T, D0, D1, D2, D3> { type Output = T; fn index(&self, ijkl:(usize,usize,usize,usize))->&T{
        let(i,j,k,l)=ijkl;
        &self.data[i][j][k][l]} 
    }
    impl<T, const D0: usize, const D1: usize, const D2: usize, const D3: usize> IndexMut<(usize,usize,usize,usize)> for Tensor4<T, D0, D1, D2, D3> { fn index_mut(&mut self, ijkl:(usize,usize,usize,usize))->&mut T{
        let(i,j,k,l)=ijkl;
        &mut self.data[i][j][k][l]
    } 
}


// -----------------------------------------------------------------------------
// Views (no allocation): lightweight slices over existing tensors
// -----------------------------------------------------------------------------
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensorView<'a, T> { 
    pub data: &'a [T], 
    pub shape: &'a [usize] 
}
#[derive(Debug, PartialEq)]
pub struct TensorViewMut<'a, T> { 
    pub data: &'a mut [T], 
    pub shape: &'a [usize] 
}


// Helper to compute flat offset for row‑major shapes
#[inline] fn offset(shape: &[usize], idx: &[usize]) -> usize {
    debug_assert_eq!(shape.len(), idx.len());
    let mut off = 0usize; let mut stride = 1usize;
    for (ax, &n) in shape.iter().rev().enumerate() {
        let i = idx[shape.len()-1-ax];
        debug_assert!(i < n);
        off += i * stride; stride *= n;
    }
    off
}

impl<'a, T> TensorView<'a, T> {
    pub fn get(&self, idx: &[usize]) -> &T { 
        &self.data[offset(self.shape, idx)] 
    }
}
impl<'a, T> TensorViewMut<'a, T> {
    pub fn get_mut(&mut self, idx: &[usize]) -> &mut T { 
        let o = offset(self.shape, idx); 
        &mut self.data[o] 
    }
}