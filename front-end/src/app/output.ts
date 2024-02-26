import { SafeUrl } from '@angular/platform-browser';

export interface RawOutput {
  text?: string
  img?: string
  alttext?: string
}

export interface Output {
  text?: string
  img?: SafeUrl
  alttext?: string;
}
