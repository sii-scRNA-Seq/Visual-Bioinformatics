import { Component, Input } from '@angular/core';
import { BlockId, BlockService } from '../block.service';

@Component({
  selector: 'app-code-block',
  templateUrl: './code-block.component.html',
  styleUrls: ['./code-block.component.css']
})
export class CodeBlockComponent {

  @Input() blockId!: BlockId;

  @Input() title: string = "";

  constructor(private blockService: BlockService) { }

  removeBlock(): void {
    this.blockService.removeBlock(this.blockId)
    console.log(this.title)
    console.log(this.blockId)
  }

}
