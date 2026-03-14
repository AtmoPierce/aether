use aether_core::math::Vector;
use aether_core::attitude::Quaternion;

#[cfg(feature = "bincode")]
use bincode::{
    config::standard,
    decode_from_slice, encode_into_slice,
    error::{DecodeError, EncodeError},
};
#[cfg(feature = "bincode")]
pub fn assert_vector_equivalence<const N: usize>(v: &Vector<f64, N>) {
    let mut buf = [0u8; 1024];    

    let len = bincode::encode_into_slice(v, &mut buf, standard())
            .expect("Buffer too small!");

    let (decoded, _): (Vector<f64, N>, usize) = bincode::decode_from_slice(&buf[..len], standard())
    .expect("Decode failed");

    assert_eq!(v, &decoded);
}

#[cfg(feature = "bincode")]
pub fn assert_quaternion_equivalence<From, To>(q: &Quaternion<f64, From, To>) 
where 
    From: aether_core::reference_frame::ReferenceFrame + std::fmt::Debug + PartialEq + bincode::Encode + bincode::Decode<()>,
    To: aether_core::reference_frame::ReferenceFrame + std::fmt::Debug + PartialEq + bincode::Encode + bincode::Decode<()>,
{
    let mut buf = [0u8; 1024];

    let len = bincode::encode_into_slice(q, &mut buf, standard())
        .expect("Failed to encode Quaternion");

    let (decoded, _): (Quaternion<f64, From, To>, usize) = 
        bincode::decode_from_slice(&buf[..len], standard())
        .expect("Failed to decode Quaternion");

    assert_eq!(q, &decoded);
}