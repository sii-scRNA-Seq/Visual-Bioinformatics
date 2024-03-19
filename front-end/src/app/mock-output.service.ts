import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable } from 'rxjs';

import { Block } from './block.interface';
import { Output } from './output';
import { OutputServiceInterface } from './output.service.interface';

@Injectable({
  providedIn: 'root'
})
export class MockOutputService implements OutputServiceInterface {
  private readonly outputs$: BehaviorSubject<Output[]> = new BehaviorSubject<Output[]> ([]);
  readonly outputs: Observable<Output[]> = this.outputs$.asObservable();
  
  async executeBlocks(blocks: Block[]): Promise<void> {
    blocks;
    const outputs = [];
    const output: Output = {text: 'Some text here'};
    outputs.push(output);
    this.outputs$.next(outputs);
  }
}
