"""
lsp -- reconstructed 1-D axial, steady-state laser-sustained argon plasma (LSP) subsystem.

Independent of the laser-delivery subsystem (the flat modules in the project root). The only
coupling is `lsp.integration`, which reads the delivery budget and builds a `LaserInput`.
"""
