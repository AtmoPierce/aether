#[derive(Debug, Clone)]
pub struct PbsTarget {
    pub architecture: Architecture,
}

#[derive(Debug, Clone, Copy)]
pub enum Architecture {
    X86_64,
    X86,
    AArch64,
    Arm,
    RiscV64,
    PowerPc64,
    S390x,
    LoongArch64,
}

impl Architecture {
    pub fn resource_value(self) -> &'static str {
        match self {
            Self::X86_64 => "x86_64",
            Self::X86 => "x86",
            Self::AArch64 => "aarch64",
            Self::Arm => "arm",
            Self::RiscV64 => "riscv64",
            Self::PowerPc64 => "powerpc64",
            Self::S390x => "s390x",
            Self::LoongArch64 => "loongarch64",
        }
    }
}

impl PbsTarget {
    pub fn from_rust_target(target: &str) -> Result<Self, String> {
        if target.contains('-') && !target.contains("-linux-") {
            return Err(format!(
                "target `{target}` is not a Linux target and cannot \
                 currently be executed by the PBS runner"
            ));
        }

        let arch = target
            .split('-')
            .next()
            .ok_or_else(|| format!("invalid Rust target: {target}"))?;

        let architecture = match arch {
            "x86_64" => Architecture::X86_64,

            "i586" | "i686" => Architecture::X86,

            "aarch64" => Architecture::AArch64,

            value if value.starts_with("arm") => Architecture::Arm,

            value if value.starts_with("riscv64") => Architecture::RiscV64,

            "powerpc64" | "powerpc64le" => Architecture::PowerPc64,

            "s390x" => Architecture::S390x,

            "loongarch64" => Architecture::LoongArch64,

            _ => {
                return Err(format!(
                    "unsupported PBS architecture for Rust target `{target}`"
                ));
            }
        };

        Ok(Self { architecture })
    }

    pub fn resource_value(&self) -> &'static str {
        self.architecture.resource_value()
    }
}

#[cfg(test)]
mod tests {
    use super::PbsTarget;

    #[test]
    fn accepts_architecture_shorthand() {
        let target = PbsTarget::from_rust_target("aarch64").unwrap();

        assert_eq!(target.resource_value(), "aarch64");
    }

    #[test]
    fn accepts_linux_rust_target() {
        let target = PbsTarget::from_rust_target("aarch64-unknown-linux-gnu").unwrap();

        assert_eq!(target.resource_value(), "aarch64");
    }

    #[test]
    fn rejects_non_linux_rust_target() {
        let error = PbsTarget::from_rust_target("aarch64-apple-darwin").unwrap_err();

        assert!(error.contains("not a Linux target"));
    }
}
