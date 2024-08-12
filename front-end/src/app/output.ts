import { SafeUrl } from '@angular/platform-browser';
import {BlockId} from './block.interface';

export interface Output {
  blockId: BlockId
  text?: string
  image?: SafeUrl
  alttext?: string
  imageList?: ImageInfo[]
}

export interface ImageInfo {
  image: SafeUrl
  alttext: string
}
