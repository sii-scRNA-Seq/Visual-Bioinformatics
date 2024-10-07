import { Observable } from 'rxjs';

import { Block } from './block.interface';
import { Output } from './output';

export interface OutputServiceInterface {

  executingBlocks: Observable<boolean>;

  outputs: Observable<Output[]>;

  resetOutputs(): void;

  executeBlocks(block: Block[]): void;
}
